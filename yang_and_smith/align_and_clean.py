#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au

"""
- Aligns the paralog fasta file using mafft, and if the option -no_supercontigs is provided, realigns using Clustal
  Omega (which can do a better job when alignment contains contigs from different regions of the full-length
  reference e.g. split between 5' and 3' halves).
- Trims alignments with Trimal.
- Runs HmmCleaner.pl on the alignments.
"""

import logging
import sys
import datetime
import textwrap
import os
import socket
import re
import glob
import subprocess
import shutil
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline, ClustalOmegaCommandline, MuscleCommandline
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait

from yang_and_smith import utils

if sys.version_info[0] < 3:
    raise Exception("Must be using Python 3")


def setup_logger(name, log_file, console_level=logging.INFO, file_level=logging.DEBUG,
                 logger_object_level=logging.DEBUG):
    """
    Function to create a logger instance.

    By default, logs level DEBUG and above to file.
    By default, logs level INFO and above to stderr and file.

    :param string name: name for the logger instance
    :param string log_file: filename for log file
    :param string console_level: level for logging to console
    :param string file_level: level for logging to file
    :param string logger_object_level: level for logger object
    :return: a logger object
    """

    # Get date and time string for log filename:
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Log to file:
    file_handler = logging.FileHandler(f'{log_file}_{date_and_time}.log', mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter('%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %('
                                    'message)s')
    file_handler.setFormatter(file_format)

    # Log to Terminal (stdout):
    console_handler = logging.StreamHandler(sys.stderr)
    console_handler.setLevel(console_level)
    console_format = logging.Formatter('%(message)s')
    console_handler.setFormatter(console_format)

    # Setup logger:
    logger_object = logging.getLogger(name)
    logger_object.setLevel(logger_object_level)  # Default level is 'WARNING'

    # Add handlers to the logger
    logger_object.addHandler(console_handler)
    logger_object.addHandler(file_handler)

    return logger_object


def mafft_or_muscle_align_multiprocessing(fasta_to_align_folder,
                                          algorithm='auto',
                                          pool_threads=1,
                                          mafft_threads=2,
                                          no_stitched_contigs=False,
                                          use_muscle=False,
                                          logger=None):
    """
    Generate alignments via function <mafft_or_muscle_align> using multiprocessing.

    :param str fasta_to_align_folder: path to folder containing input fasta files
    :param str algorithm: algorithm to use for mafft alignment; default is 'auto'
    :param int pool_threads: number of alignments to run concurrently
    :param int mafft_threads: number of threads to use for each concurrent alignment
    :param bool no_stitched_contigs: if True, realign with Clustal Omega
    :param bool use_muscle: if True, use muscle instead of mafft for alignments
    :param logging.Logger logger: a logger object
    :return str output_folder: name of the output folder containing alignments
    """

    input_folder_basename = os.path.basename(fasta_to_align_folder)
    output_folder = f'{input_folder_basename}_alignments'
    utils.createfolder(output_folder)

    if use_muscle:
        logger.info(f'{"[INFO]:":10} Generating alignments for fasta files using MUSCLE...')
    else:
        logger.info(f'{"[INFO]:":10} Generating alignments for fasta files using MAFFT...')

    # Filter out any inpout files with fewer than four sequences:
    target_genes = []
    for fasta_file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta')):
        with open(fasta_file, 'r') as input_fasta_handle:
            seqs = list(SeqIO.parse(input_fasta_handle, 'fasta'))
            if len(seqs) < 4:
                logger.info(f'{"[INFO]:":10} Skipping file {fasta_file} as it contains fewer than 4 sequences!')
                continue
            else:
                target_genes.append(fasta_file)

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(mafft_or_muscle_align,
                                      fasta_file,
                                      algorithm,
                                      output_folder,
                                      counter, lock,
                                      num_files_to_process=len(target_genes),
                                      threads=mafft_threads,
                                      no_stitched_contigs=no_stitched_contigs,
                                      use_muscle=use_muscle,
                                      logger=logger)

                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(utils.done_callback)
        wait(future_results, return_when="ALL_COMPLETED")

    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                      utils.file_exists_and_not_empty(alignment)]

    logger.debug(f'{"[INFO]:":10} {len(alignment_list)} alignments generated from {len(future_results)} fasta files...')

    return output_folder


def mafft_or_muscle_align(fasta_file,
                          algorithm,
                          output_folder,
                          counter,
                          lock,
                          num_files_to_process,
                          threads=1,
                          no_stitched_contigs=False,
                          use_muscle=False,
                          logger=None):
    """
    Use mafft or muscle to align a fasta file of sequences, using the algorithm (if mafft) and number of threads
    provided. Trims alignment with Trimal if no_stitched_contigs=False. Returns filename of the output alignment.

    :param str fasta_file: path to a fasta file
    :param str algorithm: algorithm to use for mafft alignment; default is 'auto'
    :param str output_folder: name of output folder for alignments
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int num_files_to_process: total number of fasta files for alignment
    :param int threads: number of threads to use for alignment program
    :param bool no_stitched_contigs: if True, realign with Clustal Omega
    :param bool use_muscle: if True, use muscle instead of mafft for alignments
    :param logging.Logger logger: a logger object
    :return str expected_alignment_file/expected_alignment_file_trimmed: filename of output alignment
    """

    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'
    expected_alignment_file_trimmed = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)

    try:
        if not no_stitched_contigs:  # i.e. stitched contig _were_ produced
            assert utils.file_exists_and_not_empty(expected_alignment_file_trimmed)
            logger.debug(f'Trimmed alignment exists for {fasta_file_basename}, skipping...')
            with lock:
                counter.value += 1

            return os.path.basename(expected_alignment_file_trimmed)
        else:
            assert utils.file_exists_and_not_empty(expected_alignment_file)
            logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
            with lock:
                counter.value += 1

            return os.path.basename(expected_alignment_file)

    except AssertionError:
        if use_muscle:
            logger.info(f'{"[INFO]:":10} Alignment will be performed using MUSCLE rather than MAFFT!')
            muscle_cline = MuscleCommandline(input=fasta_file, out=expected_alignment_file)
            stdout, stderr = muscle_cline()

            logger.debug(f'stdout is: {stdout}')
            logger.debug(f'stderr is: {stderr}')
        else:
            if algorithm == 'auto':
                mafft_cline = (MafftCommandline(auto='true', adjustdirection='true', thread=threads, input=fasta_file))
            else:
                mafft_cline = (MafftCommandline(algorithm, adjustdirection='true', thread=threads, input=fasta_file))

            logger.info(f'{"[INFO]:":10} Performing MAFFT alignment with command: {mafft_cline}')
            stdout, stderr = mafft_cline()

            logger.debug(f'stdout is: {stdout}')
            logger.debug(f'stderr is: {stderr}')

            with open(expected_alignment_file, 'w') as alignment_file:
                alignment_file.write(stdout)

            # If mafft has reversed any sequences, remove the prefix "_R_" from such sequence names:
            remove_r_prefix(expected_alignment_file, logger=logger)

        if not no_stitched_contigs:  # only trim if stitched contigs (and hence won't be realigned)
            trimmed_alignment = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)
            try:
                result = subprocess.run(['trimal', '-in', expected_alignment_file, '-out', trimmed_alignment,
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
                sys.exit('There was an issue running trimal. Check input files!')

        with lock:
            counter.value += 1
            logger.debug(f'Aligned file {fasta_file_basename}')

        return os.path.basename(expected_alignment_file)

    finally:
        with lock:
            sys.stderr.write(f'\r{"[INFO]:":10} Finished generating alignment for file {fasta_file_basename},'
                             f' {counter.value}/{num_files_to_process}')


def remove_r_prefix(alignment, logger=None):
    """
    Takes a fasta alignment, removes any '_R_' prefix in fasta headers (inserted by mafft if a sequences was
    reversed) and writes a new alignment to the same filename (i.e overwrites the original file).

    :param str alignment: oath to untrimmed alignment fasta file
    :param logging.Logger logger: a logger object
    :return:
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
        sequence_names = ', '.join(seqs_renamed)

        fill = textwrap.fill(f'{"[WARNING]:":10} Alignment {alignment} contains sequences that were reversed by '
                             f'MAFFT, i.e. sequence names have the prefix "_R_". This prefix has been removed from the '
                             f'following sequences: {sequence_names}',
                             width=90,
                             subsequent_indent=" " * 11)
        logger.warning(fill)


def run_hmm_cleaner(input_folder, logger=None):
    """
    Runs HmmCleaner.pl on each alignment within a provided folder.

    :param str input_folder: path to a folder containing trimmed fasta alignment files
    :param logging.Logger logger: a logger object
    :return:
    """

    input_folder_basename = os.path.basename(input_folder)
    output_folder = f'{input_folder_basename}_hmmcleaned'
    utils.createfolder(output_folder)

    logger.info(f'\n{"[INFO]:":10} Running HmmCleaner.pl on trimmed alignments...')

    for alignment in glob.glob(f'{input_folder}/*.aln.trimmed.fasta'):
        command = f'/usr/bin/perl /usr/local/bin/HmmCleaner.pl {alignment}'

        host = socket.gethostname()
        if host == 'RBGs-MacBook-Air.local':
            command = f'/Users/chrisjackson/perl5/perlbrew/perls/perl-5.26.2/bin/perl ' \
                      f'/Users/chrisjackson/perl5/perlbrew/perls/perl-5.26.2/bin/HmmCleaner.pl {alignment}'

        logger.debug(f'trying command {command}')

        try:
            result = subprocess.run(command, shell=True, universal_newlines=True, check=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)
            logger.debug(f'hmmcleaner check_returncode() is: {result.check_returncode()}')
            logger.debug(f'hmmcleaner stdout is: {result.stdout}')
            logger.debug(f'hmmcleaner stderr is: {result.stderr}')

            # Filter out empty sequences comprised only of dashes, and post-hmmcleaner alignments where all sequences
            # are either dashes or empty. If fewer than 4 'good' sequences are present, skip the gene:
            hmm_file = re.sub('aln.trimmed.fasta', 'aln.trimmed_hmm.fasta', str(alignment))
            hmm_file_output = re.sub('aln.trimmed.fasta', 'aln.hmm.trimmed.fasta', str(alignment))
            with open(hmm_file, 'r') as hmm_fasta:
                good_seqs = []
                seqs_all_dashes = []
                empty_seqs = []
                seqs = SeqIO.parse(hmm_fasta, 'fasta')
                for seq in seqs:
                    characters = set(character for character in seq.seq)
                    if len(characters) == 0:
                        empty_seqs.append(seq.name)
                    elif len(characters) == 1 and '-' in characters:
                        seqs_all_dashes.append(seq.name)
                    else:
                        good_seqs.append(seq)

            if seqs_all_dashes:
                seqs_all_dashes_joined = ', '.join(seqs_all_dashes)
                logger.debug(f'After running HmmCleaner.pl, the following sequences contained dashes, and have been '
                             f'removed: {seqs_all_dashes_joined}')

            if empty_seqs:
                empty_seqs_joined = ', '.join(empty_seqs)
                logger.debug(f'After running HmmCleaner.pl, the following sequences were empty, and have been '
                             f'removed: {empty_seqs_joined}')

            if len(good_seqs) < 4:
                logger.warning(f'{"[WARNING]:":10} After running HmmCleaner.pl, file {os.path.basename(hmm_file)} '
                               f'contains fewer than 4 good sequences, skipping gene!')
            else:
                with open(hmm_file_output, 'w') as filtered_hmm_fasta:
                    SeqIO.write(good_seqs, filtered_hmm_fasta, 'fasta')

        except subprocess.CalledProcessError as exc:
            logger.error(f'hmmcleaner FAILED. Output is: {exc}')
            logger.error(f'hmmcleaner stdout is: {exc.stdout}')
            logger.error(f'hmmcleaner stderr is: {exc.stderr}')

            hmm_file_output = re.sub('aln.trimmed.fasta', 'aln.hmm.trimmed.fasta', str(alignment))

            logger.info(f'{"[INFO]:":10} Could not run HmmCleaner.pl for alignment {alignment} using command'
                        f' {command}')
            logger.info(f'Copying alignment {alignment} to {hmm_file_output} anyway...')

            shutil.copy(alignment, hmm_file_output)

    # Copy post-hmmcleaner (successful or not) files to a new output directory:
    for file in glob.glob(f"{input_folder}/*aln.hmm.trimmed*"):
        try:
            shutil.move(file, output_folder)
        except shutil.Error as e:
            logger.debug(f'Copying file {file} produced the error: {e}')


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

    logger.debug(f'{"[INFO]:":10} {len(alignment_list)} alignments generated from {len(future_results)} fasta files...')

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

        # remove_r_prefix(expected_alignment_file)
        trimmed_alignment = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)

        try:
            result = subprocess.run(['trimal', '-in', expected_alignment_file, '-out', trimmed_alignment,
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
            sys.exit('There was an issue running trimal. Check input files!')

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

def main(args):
    """
    Entry point for the resolve_paralogs.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = setup_logger(__name__, 'align_and_clean')

    logger.info(f'{"[INFO]:":10} Script was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')

    if not args.no_stitched_contigs:  # i.e. if it's a standard run with stitched contigs produced.
        logger.debug(f'Running without no_stitched_contigs option - aligning with mafft or muscle only')
        alignments_output_folder = mafft_or_muscle_align_multiprocessing(
            args.gene_fasta_directory,
            algorithm=args.mafft_algorithm,
            pool_threads=args.pool,
            mafft_threads=args.threads,
            no_stitched_contigs=args.no_stitched_contigs,
            use_muscle=args.use_muscle,
            logger=logger)

        run_hmm_cleaner(alignments_output_folder, logger=logger)

    elif args.no_stitched_contigs:  # Re-align with Clustal Omega.
        logger.debug(f'Running with no_stitched_contigs option - realigning with clustal omega')
        alignments_output_folder = mafft_or_muscle_align_multiprocessing(
            args.gene_fasta_directory,
            algorithm=args.mafft_algorithm,
            pool_threads=args.pool,
            mafft_threads=args.threads,
            no_stitched_contigs=args.no_stitched_contigs,
            use_muscle=args.use_muscle,
            logger=logger)

        clustal_alignment_output_folder = clustalo_align_multiprocessing(
            alignments_output_folder,
            pool_threads=args.pool,
            clustalo_threads=args.threads,
            logger=logger)

        run_hmm_cleaner(clustal_alignment_output_folder, logger=logger)

