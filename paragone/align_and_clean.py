#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
- Aligns the paralog fasta files using MAFFT, and if the option -no_stitched_contigs is provided,
  realigns using Clustal Omega (which can do a better job when the alignment contains contigs from different regions of
  the full-length reference e.g. split between 5' and 3' halves).
- Trims alignments with Trimal (optional)
- Runs HmmCleaner.pl on the alignments (optional).
"""

import logging
import sys
import textwrap
import os
import socket
import re
import glob
import subprocess
import shutil
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline, ClustalOmegaCommandline
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait

from paragone import utils


def mafft_align_multiprocessing(fasta_to_align_folder,
                                algorithm='auto',
                                pool_threads=1,
                                mafft_threads=2,
                                logger=None):
    """
    Generate alignments via function <mafft_align> using multiprocessing.

    :param str fasta_to_align_folder: path to folder containing input fasta files
    :param str algorithm: algorithm to use for mafft alignment; default is 'auto'
    :param int pool_threads: number of alignments to run concurrently
    :param int mafft_threads: number of threads to use for each concurrent alignment
    :param logging.Logger logger: a logger object
    :return str output_folder: name of the output folder containing alignments
    """

    output_folder = f'02_alignments'
    utils.createfolder(output_folder)

    logger.info(f'{"[INFO]:":10} Generating alignments for fasta files using MAFFT...')

    # Filter out any input files with fewer than four sequences:
    target_genes = []
    for fasta_file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta')):
        with open(fasta_file, 'r') as input_fasta_handle:
            seqs = list(SeqIO.parse(input_fasta_handle, 'fasta'))
            if len(seqs) < 4:
                logger.warning(f'{"[WARNING]:":10} Skipping file {fasta_file} as it contains fewer than 4 sequences!')
                continue
            else:
                target_genes.append(fasta_file)

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(mafft_align,
                                      fasta_file,
                                      algorithm,
                                      output_folder,
                                      counter,
                                      lock,
                                      num_files_to_process=len(target_genes),
                                      threads=mafft_threads,
                                      logger=logger)

                          for fasta_file in target_genes]

        for future in future_results:
            future.add_done_callback(utils.done_callback)
        wait(future_results, return_when="ALL_COMPLETED")

    for future in future_results:
        try:
            alignment_name, seqs_renamed = future.result()
            if seqs_renamed:
                sequence_names = ', '.join(seqs_renamed)
                logger.info('')
                fill = textwrap.fill(f'{"[WARNING]:":10} Alignment {alignment_name} contains sequences that were '
                                     f'reversed by MAFFT, i.e. sequence names have the prefix "_R_". This prefix has '
                                     f'been removed from the following sequences: {sequence_names}',
                                     width=90, subsequent_indent=" " * 11)
                logger.warning(fill)
        except:
            raise

    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                      utils.file_exists_and_not_empty(alignment)]

    logger.debug(f'{len(alignment_list)} alignments generated from {len(future_results)} fasta files...')

    return output_folder


def mafft_align(fasta_file,
                algorithm,
                output_folder,
                counter,
                lock,
                num_files_to_process,
                threads=1,
                logger=None):
    """
    Use mafft to align a fasta file of sequences, using the algorithm (default is 'auto') and number of threads
    provided. Returns filename of the output alignment.

    :param str fasta_file: path to a fasta file
    :param str algorithm: algorithm to use for mafft alignment; default is 'auto'
    :param str output_folder: name of output folder for alignments
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int num_files_to_process: total number of fasta files for alignment
    :param int threads: number of threads to use for alignment program
    :param logging.Logger logger: a logger object
    :return str expected_alignment_file/expected_alignment_file_trimmed: filename of output alignment
    """

    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'

    try:
        assert utils.file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
        with lock:
            counter.value += 1

        return os.path.basename(expected_alignment_file), None  # None here as it should already have been processed

    except AssertionError:

        if algorithm == 'auto':
            mafft_cline = (MafftCommandline(auto='true', adjustdirection='true', thread=threads, input=fasta_file))
        else:
            mafft_cline = (MafftCommandline(algorithm, adjustdirection='true', thread=threads, input=fasta_file))

        logger.debug(f'{"[INFO]:":10} Performing MAFFT alignment with command: {mafft_cline}')

        stdout, stderr = mafft_cline()
        logger.debug(f'stdout is: {stdout}')
        logger.debug(f'stderr is: {stderr}')

        with open(expected_alignment_file, 'w') as alignment_file:
            alignment_file.write(stdout)

        # If mafft has reversed any sequences, remove the prefix "_R_" from such sequence names:
        seqs_renamed_list = remove_r_prefix(expected_alignment_file, logger=logger)

        with lock:
            counter.value += 1
            logger.debug(f'Aligned file {fasta_file_basename}')

        return os.path.basename(expected_alignment_file), seqs_renamed_list

    finally:
        with lock:
            sys.stderr.write(f'\r{"[INFO]:":10} Finished generating alignment for file {fasta_file_basename},'
                             f' {counter.value}/{num_files_to_process}')


def remove_r_prefix(alignment, logger=None):
    """
    Takes a fasta alignment, removes any '_R_' prefix in fasta headers (inserted by mafft if a sequence was
    reversed) and writes a new alignment to the same filename (i.e overwrites the original file). Note that it is
    expected that mafft might reverse some sequences for correct alignment, as the paralog fasta files output by
    HybPiper can contain revcomp sequences.

    :param str alignment: path to untrimmed alignment fasta file
    :param logging.Logger logger: a logger object
    :return str os.path.basename(alignment), list seqs_renamed:
    """

    seqs_renamed = []
    with open(alignment) as alignment_handle:
        alignment_obj = AlignIO.read(alignment_handle, 'fasta')
        for seq in alignment_obj:
            if seq.name.startswith('_R_'):
                seqs_renamed.append(seq.name)
                seq.name = seq.name.lstrip('_R_')
                seq.id = seq.id.lstrip('_R_')
        with open(alignment, 'w') as new_alignment_handle:
            AlignIO.write(alignment_obj, new_alignment_handle, 'fasta')

    if seqs_renamed:
        return seqs_renamed
    else:
        return None


def run_trimal(input_folder,
               output_folder,
               logger=None):
    """
    Runs trimal on each alignment within a provided folder.
    :param str input_folder: path to a folder containing fasta alignment files
    :param str output_folder: path to a folder for output trimmed fasta alignment files
    :param logging.Logger logger: a logger object
    :return str trimmed_alignments_directory: directory containing trimmed alignments
    """

    trimmed_alignments_directory = utils.createfolder(output_folder)

    logger.info('')
    fill = textwrap.fill(f'{"[INFO]:":10} Running trimal on paralog alignments. Trimmed alignments will be '
                         f'written to directory: "{output_folder}".',
                         width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
    logger.info(fill)

    for alignment in glob.glob(f'{input_folder}/*.aln.fasta'):
        alignment_basename = os.path.basename(alignment)
        expected_trimmed_alignment_file = \
            f'{output_folder}/{re.sub(".aln.fasta", ".aln.trimmed.fasta", alignment_basename)}'

        try:
            assert utils.file_exists_and_not_empty(expected_trimmed_alignment_file)
            logger.debug(f'Trimmed alignment exists for {alignment_basename}, skipping...')

            continue

        except AssertionError:
            try:
                result = subprocess.run(['trimal', '-in', alignment, '-out', expected_trimmed_alignment_file,
                                         '-gapthreshold', '0.12', '-terminalonly', '-gw', '1'],
                                        universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        check=True)
                logger.debug(f'trimal check_returncode() is: {result.check_returncode()}')
                logger.debug(f'trimal stdout is: {result.stdout}')
                logger.debug(f'trimal stderr is: {result.stderr}')

            except subprocess.CalledProcessError as exc:
                logger.error(f'trimal FAILED. Output is: {exc}')
                logger.error(f'trimal stdout is: {exc.stdout}')
                logger.error(f'trimal stderr is: {exc.stderr}')
                raise ValueError('There was an issue running trimal. Check input files!')

    return trimmed_alignments_directory


def run_hmm_cleaner(input_folder,
                    no_trimming=False,
                    logger=None):
    """
    Runs HmmCleaner.pl on each alignment within a provided folder.

    :param str input_folder: path to a folder containing trimmed fasta alignment files
    :param bool no_trimming: if True, sequences have NOT been trimmed with trimal
    :param logging.Logger logger: a logger object
    :return:
    """

    if no_trimming:
        output_folder = f'04_alignments_hmmcleaned'
    else:
        output_folder = f'04_alignments_trimmed_hmmcleaned'

    utils.createfolder(output_folder)

    logger.info('')
    fill = textwrap.fill(f'{"[INFO]:":10} Running HmmCleaner.pl on alignments. Cleaned alignments will be '
                         f'written to directory: "{output_folder}".',
                         width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
    logger.info(fill)

    for alignment in glob.glob(f'{input_folder}/*.fasta'):

        alignment_basename = os.path.basename(alignment)
        command = f'perl $(which HmmCleaner.pl) {alignment}'

        logger.debug(f'Trying command {command}')

        hmm_file = re.sub('.fasta', '_hmm.fasta', str(alignment_basename))  # output by HmmCleaner.pl
        hmm_score = re.sub('.fasta', '_hmm.score', str(alignment_basename))  # output by HmmCleaner.pl
        hmm_log = re.sub('.fasta', '_hmm.log', str(alignment_basename))  # output by HmmCleaner.pl
        hmm_file_output = re.sub('.fasta', '.hmm.fasta', str(alignment_basename))  # Desired filename

        try:
            result = subprocess.run(command, shell=True, universal_newlines=True, check=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            logger.debug(f'hmmcleaner check_returncode() is: {result.check_returncode()}')
            logger.debug(f'hmmcleaner stdout is: {result.stdout}')
            logger.debug(f'hmmcleaner stderr is: {result.stderr}')

            # Filter out empty sequences comprising only dashes, and post-HmmCleaner alignments where all sequences
            # are either dashes or empty. If fewer than 4 'good' sequences are present, skip the gene:
            with open(f'{input_folder}/{hmm_file}', 'r') as hmm_fasta_handle:
                good_seqs = []
                seqs_all_dashes = []
                empty_seqs = []
                seqs = SeqIO.parse(hmm_fasta_handle, 'fasta')

                for seq in seqs:
                    characters = set(character for character in seq.seq)
                    if len(characters) == 0:
                        empty_seqs.append(seq.name)
                    elif len(characters) == 1 and '-' in characters:
                        seqs_all_dashes.append(seq.name)
                    else:
                        good_seqs.append(seq)

            # Log any sequences that were removed:
            if seqs_all_dashes:
                seqs_all_dashes_joined = ', '.join(seqs_all_dashes)
                logger.debug(f'After running HmmCleaner.pl, the following sequences contained dashes, and have been '
                             f'removed: {seqs_all_dashes_joined}')

            if empty_seqs:
                empty_seqs_joined = ', '.join(empty_seqs)
                logger.debug(f'After running HmmCleaner.pl, the following sequences were empty, and have been '
                             f'removed: {empty_seqs_joined}')

            # Skip any filtered alignments with fewer than 4 sequences remaining:
            if len(good_seqs) < 4:
                logger.warning(f'{"[WARNING]:":10} After running HmmCleaner.pl, file {os.path.basename(hmm_file)} '
                               f'contains fewer than 4 good sequences, skipping gene!')

            else:
                with open(f'{output_folder}/{hmm_file_output}', 'w') as filtered_hmm_fasta:
                    SeqIO.write(good_seqs, filtered_hmm_fasta, 'fasta')

                # Remove the original HmmCleaner output fasta file:
                os.remove(f'{input_folder}/{hmm_file}')

                # Delete the HmmCleaner score and log files:
                os.remove(f'{input_folder}/{hmm_score}')
                os.remove(f'{input_folder}/{hmm_log}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'hmmcleaner FAILED. Output is: {exc}')
            logger.error(f'hmmcleaner stdout is: {exc.stdout}')
            logger.error(f'hmmcleaner stderr is: {exc.stderr}')

            logger.info(f'{"[INFO]:":10} Could not run HmmCleaner.pl for alignment {alignment} using command'
                        f' {command}')
            logger.info(f'Copying un-cleaned alignment {alignment} to {hmm_file_output} anyway...')

            shutil.copy(alignment, f'{output_folder}/{hmm_file_output}')

        except:
            raise


def clustalo_align_multiprocessing(fasta_to_align_folder,
                                   pool_threads=1,
                                   clustalo_threads=1,
                                   logger=None):
    """
    Generate alignments via function <clustalo_align> using multiprocessing.

    :param fasta_to_align_folder: path to folder containing input fasta alignment files
    :param int pool_threads: number of alignments to run concurrently
    :param int clustalo_threads: number of threads to use for each concurrent alignment
    :param logging.Logger logger: a logger object
    :return str output_folder: name of the output folder containing alignments
    """

    input_folder_basename = os.path.basename(fasta_to_align_folder)
    output_folder = f'{input_folder_basename}_clustal'
    utils.createfolder(output_folder)

    logger.info(f'{"[INFO]:":10} Generating alignments for fasta files using Clustal Omega...')
    target_genes = [file for file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(clustalo_align,
                                      fasta_file,
                                      output_folder,
                                      counter,
                                      lock,
                                      num_files_to_process=len(target_genes),
                                      threads=clustalo_threads,
                                      logger=logger)
                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(utils.done_callback)
        wait(future_results, return_when="ALL_COMPLETED")

    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                      utils.file_exists_and_not_empty(alignment)]

    logger.debug(f'{len(alignment_list)} alignments generated from {len(future_results)} fasta files...')

    return output_folder


def clustalo_align(fasta_file,
                   output_folder,
                   counter,
                   lock,
                   num_files_to_process,
                   threads=1,
                   logger=None):
    """
    Use clustal omega to align a fasta file of sequences, using the number of threads provided. Trims alignment with
    Trimal. Returns filename of the alignment produced.

    :param str fasta_file: path to a fasta file
    :param str output_folder: name of output folder for alignments
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int num_files_to_process: total number of fasta files for alignment
    :param int threads: number of threads to use for alignment program
    :param logging.Logger logger: a logger object
    :return str expected_alignment_file: filename of output alignment
    """

    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{fasta_file_basename}'

    try:
        assert utils.file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
        with lock:
            counter.value += 1

        return os.path.basename(expected_alignment_file)

    except AssertionError:
        clustalomega_cline = ClustalOmegaCommandline(infile=fasta_file, outfile=expected_alignment_file,
                                                     verbose=True, auto=True, threads=threads)

        logger.info(f'{"[INFO]:":10} Performing Clustal alignment with command: {clustalomega_cline}')

        stdout, stderr = clustalomega_cline()

        logger.debug(f'stdout is: {stdout}')
        logger.debug(f'stderr is: {stderr}')

        with lock:
            counter.value += 1
            logger.debug(f'Aligned file {fasta_file_basename}')

        return os.path.basename(expected_alignment_file)

    finally:
        with lock:
            sys.stderr.write(f'\r{"[INFO]:":10} Finished generating alignment for file {fasta_file_basename},'
                             f' {counter.value}/{num_files_to_process}')


########################################################################################################################
########################################################################################################################
# Run script:

def main(args, logger=None):
    """
    Entry point for the paragone_main.py script

    :param args: argparse namespace with subparser options for function main()
    :param logging.Logger logger: a logger object
    :return:
    """

    logger.debug(f'{"[INFO]:":10} Module align_and_clean was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)

    logger.debug(f'{fill}')
    logger.debug(args)

    logger.info(f'{"[INFO]:":10} ======> ALIGNING PARALOGS AND TRIMMING/CLEANING ALIGNMENTS <======\n')

    gene_fasta_directory = '01_input_paralog_fasta_with_sanitised_filenames'

    # Checking input directories and files:
    directory_suffix_dict = {gene_fasta_directory: '.fasta'}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    if not args.no_stitched_contigs:  # i.e. if it's a standard run with stitched contigs produced.
        logger.debug(f'Running without no_stitched_contigs option - aligning with mafft only')

        alignments_output_folder = mafft_align_multiprocessing(
            gene_fasta_directory,
            algorithm=args.mafft_algorithm,
            pool_threads=args.pool,
            mafft_threads=args.threads,
            logger=logger)

        # Perform optional trimming with TrimAl:
        trimmed_output_folder = '03_alignments_trimmed'
        if not args.no_trimming:
            alignments_output_folder = run_trimal(alignments_output_folder,
                                                  trimmed_output_folder,
                                                  logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping trimming step...')

        # Perform optional cleaning with HmmCleaner.pl
        if not args.no_cleaning:
            run_hmm_cleaner(alignments_output_folder,
                            no_trimming=args.no_trimming,
                            logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping cleaning step...')

    elif args.no_stitched_contigs:  # Re-align with Clustal Omega.
        logger.debug(f'Running with no_stitched_contigs option - realigning with clustal omega')

        # Align with MAFFT first to reverse complement any paralog sequences that need it:
        alignments_output_folder = mafft_align_multiprocessing(
            gene_fasta_directory,
            algorithm=args.mafft_algorithm,
            pool_threads=args.pool,
            mafft_threads=args.threads,
            logger=logger)

        alignments_output_folder = clustalo_align_multiprocessing(
            alignments_output_folder,
            pool_threads=args.pool,
            clustalo_threads=args.threads,
            logger=logger)

        # Perform optional trimming with TrimAl:
        trimmed_output_folder = '03_alignments_trimmed'
        if not args.no_trimming:
            alignments_output_folder = run_trimal(alignments_output_folder,
                                                  trimmed_output_folder,
                                                  logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping trimming step...')

        # Perform optional cleaning with HmmCleaner.pl
        if not args.no_cleaning:
            run_hmm_cleaner(alignments_output_folder,
                            no_trimming=args.no_trimming,
                            logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping cleaning step...')

    logger.info(f'{"[INFO]:":10} Finished aligning and trimming/cleaning input files.')

