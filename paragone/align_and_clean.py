#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
- Aligns the paralog fasta files using MAFFT, and if the option -no_stitched_contigs is provided,
  realigns using Clustal Omega (which can do a better job when the alignment contains contigs from different regions of
  the full-length reference e.g. split between 5' and 3' halves).
- Trims alignments with Trimal (optional)
- Runs TAPER (https://github.com/chaoszhang/TAPER on the alignments (optional).
"""

import logging
import sys
import textwrap
import os
import re
import glob
import subprocess
import io
import traceback
from Bio import SeqIO, AlignIO
from concurrent.futures.process import ProcessPoolExecutor
from concurrent.futures import as_completed
from multiprocessing import Manager
from concurrent.futures import wait

from paragone import utils


def mafft_align_multiprocessing(fasta_to_align_folder,
                                algorithm='auto',
                                adjust_direction=False,
                                pool_threads=1,
                                mafft_threads=2,
                                logger=None):
    """
    Generate alignments via function <mafft_align> using multiprocessing.

    https://mafft.cbrc.jp/alignment/software/adjustdirection.html

    :param str fasta_to_align_folder: path to folder containing input fasta files
    :param str algorithm: algorithm to use for mafft alignment; default is 'auto'
    :param bool adjust_direction: if True, enable the --adjustdirection flag in MAFFT
    :param int pool_threads: number of alignments to run concurrently
    :param int mafft_threads: number of threads to use for each concurrent alignment
    :param logging.Logger logger: a logger object
    :return str output_folder: name of the output folder containing alignments
    """

    output_folder = f'02_alignments'
    utils.createfolder(output_folder)

    logger.info(f'{"[INFO]:":10} Generating alignments for fasta files using MAFFT...')

    if adjust_direction:
        fill = textwrap.fill(f'{"[INFO]:":10} The MAFFT flag "--adjustdirection" is set. This allows MAFFT to '
                             f'generate reverse complement sequences, as necessary, and align them together with the '
                             f'remaining sequences. Note that MAFFT assumes that the first sequence in the input '
                             f'*.fasta file is in the correct orientation!',
                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

        logger.info(f'{fill}')

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
                                      adjust_direction,
                                      output_folder,
                                      counter,
                                      lock,
                                      num_files_to_process=len(target_genes),
                                      threads=mafft_threads,
                                      logger=logger)

                          for fasta_file in target_genes]

        for future in as_completed(future_results):

            try:
                check = future.result()

            except Exception as error:
                logger.error(f'Error raised: {error}')
                tb = traceback.format_exc()
                logger.error(f'traceback is:\n{tb}')
                sys.exit(1)

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
                adjust_direction,
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
    :param bool adjust_direction: if True, enable the --adjustdirection flag in MAFFT
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

        return os.path.basename(expected_alignment_file), None  # "None" here as it should already have been processed

    except AssertionError:
        if algorithm == 'auto':
            if adjust_direction:
                command = f'mafft --auto --thread {threads} --adjustdirection {fasta_file} > {expected_alignment_file}'
            else:
                command = f'mafft --auto --thread {threads} {fasta_file} > {expected_alignment_file}'
        else:
            if adjust_direction:
                command = f'{algorithm} --thread {threads} --adjustdirection {fasta_file} > {expected_alignment_file}'
            else:
                command = f'{algorithm} --thread {threads} {fasta_file} > {expected_alignment_file}'

        logger.debug(f'{"[INFO]:":10} Performing MAFFT alignment with command: {command}')

        try:

            result = subprocess.run(command,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    check=True,
                                    shell=True)

            logger.debug(f'MAFFT check_returncode() is: {result.check_returncode()}')
            logger.debug(f'MAFFT stdout is: {result.stdout}')
            logger.debug(f'MAFFT stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'MAFFT FAILED. Output is: {exc}')
            logger.error(f'MAFFT stdout is: {exc.stdout}')
            logger.error(f'MAFFT stderr is: {exc.stderr}')
            raise ValueError('There was an issue running MAFFT. Check input files!')

        # If mafft has reversed any sequences, remove the prefix "_R_" from corresponding sequence names:
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
               args_for_trimal_options=None,
               logger=None):
    """
    Runs trimal on each alignment within a provided folder.
    :param str input_folder: path to a folder containing fasta alignment files
    :param str output_folder: path to a folder for output trimmed fasta alignment files
    :param argparse.Namespace args_for_trimal_options: argparse.Namespace object to get trimal options, else None
    :param logging.Logger logger: a logger object
    :return str trimmed_alignments_directory: directory containing trimmed alignments
    """

    trimal_options_string = utils.get_trimal_options(args_for_trimal_options,
                                                     logger=logger)

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
                command = (f'trimal '
                           f'-in {alignment} '
                           f' {trimal_options_string}')

                logger.debug(f'trimal command is: {command}')

                result = subprocess.run(command,
                                        universal_newlines=True,
                                        stdout=subprocess.PIPE,
                                        stderr=subprocess.PIPE,
                                        check=True,
                                        shell=True)

                logger.debug(f'trimal check_returncode() is: {result.check_returncode()}')
                logger.debug(f'trimal stdout is: {result.stdout}')
                logger.debug(f'trimal stderr is: {result.stderr}')

            except subprocess.CalledProcessError as exc:
                logger.error(f'trimal FAILED. Output is: {exc}')
                logger.error(f'trimal stdout is: {exc.stdout}')
                logger.error(f'trimal stderr is: {exc.stderr}')
                raise ValueError('There was an issue running trimal. Check input files!')

        # Filter out empty sequences comprising only dashes, and post-Trimal alignments where all sequences are either
        # dashes or empty. If fewer than 4 'good' sequences are present, skip the alignment:
        good_seqs = []
        seqs_all_dashes = []
        empty_seqs = []

        alignment_io = io.StringIO(result.stdout)
        seqs = SeqIO.parse(alignment_io, 'fasta')

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

            fill = textwrap.fill(f'{"[INFO]:":10} After running Trimal, the following sequences in alignment'
                                 f' {alignment_basename} contained only dashes, and have been removed:',
                                 width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
            logger.info(fill)
            logger.info(f'{" "* 11}{seqs_all_dashes_joined}')

        if empty_seqs:
            empty_seqs_joined = ', '.join(empty_seqs)

            fill = textwrap.fill(f'{"[INFO]:":10} After running Trimal, the following sequences in alignment'
                                 f' {alignment_basename} contained were empty, and have been removed:',
                                 width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
            logger.info(fill)
            logger.info(f'{" " * 11}{empty_seqs_joined}')

        # Skip any filtered alignments with fewer than 4 sequences remaining:
        if len(good_seqs) < 4:

            fill = textwrap.fill(f'{"[INFO]:":10} After running Trimal, alignment'
                                 f' {alignment_basename} contains fewer than 4 sequences, skipping alignment!',
                                 width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
            logger.info(fill)

        else:

            with open(f'{expected_trimmed_alignment_file}', 'w') as trimal_checked_fasta:
                SeqIO.write(good_seqs, trimal_checked_fasta, 'fasta')

    return trimmed_alignments_directory


def run_taper_multiprocessing(alignments_to_clean_folder,
                              no_trimming=False,
                              cleaning_cutoff=3,
                              pool_threads=1,
                              logger=None):
    """
    Runs TAPER (correction_multi.jl) on each alignment within a provided folder.

    :param str alignments_to_clean_folder: path to folder containing input fasta alignment files for cleaning
    :param bool no_trimming: if True, sequences have NOT been trimmed with trimal
    :param int cleaning_cutoff: cutoff value to pass to TAPER. Lower will perform more aggressive cleaning
    :param int pool_threads: number of alignments to clean concurrently
    :param logging.Logger logger: a logger object
    :return str output_folder: name of the output folder containing cleaned alignments
    """

    if no_trimming:
        output_folder = f'04_alignments_cleaned'
    else:
        output_folder = f'04_alignments_trimmed_cleaned'

    utils.createfolder(output_folder)

    logger.info('')
    fill = textwrap.fill(f'{"[INFO]:":10} Running TAPER on alignments. Cleaned alignments will be '
                         f'written to directory: "{output_folder}".',
                         width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
    logger.info(fill)

    input_alignments = [file for file in sorted(glob.glob(f'{alignments_to_clean_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(run_taper,
                                      fasta_file,
                                      output_folder,
                                      counter,
                                      lock,
                                      cleaning_cutoff=cleaning_cutoff,
                                      num_files_to_process=len(input_alignments),
                                      logger=logger)
                          for fasta_file in input_alignments]

        for future in as_completed(future_results):

            try:
                check = future.result()

                # Log any messages from the TAPER processes to the main log:
                cleaned_alignment, log_list_check_seqs = check
                if log_list_check_seqs:
                    for item in log_list_check_seqs:
                        fill = textwrap.fill(f'{"[INFO]:":10} {item}',
                                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
                        logger.info(fill)

            except Exception as error:
                logger.error(f'Error raised: {error}')
                tb = traceback.format_exc()
                logger.error(f'traceback is:\n{tb}')
                sys.exit(1)

        wait(future_results, return_when="ALL_COMPLETED")  # redundant, but...

    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.hmm.fasta') if
                      utils.file_exists_and_not_empty(alignment)]

    logger.debug(f'{len(alignment_list)} Cleaned alignments generated from {len(future_results)} alignment files...')

    return output_folder


def run_taper(alignment,
              output_folder,
              counter,
              lock,
              cleaning_cutoff,
              num_files_to_process,
              logger=None):
    """
    Runs TAPER on an alignment.

    :param str alignment: name of a fasta alignment file
    :param str output_folder: name of output folder for cleaned sequences
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int cleaning_cutoff: cutoff value to pass to TAPER. Lower will perform more aggressive cleaning
    :param int num_files_to_process: total number of fasta files for alignment
    :param logging.Logger logger: a logger object
    :return:
    """

    alignment_basename = os.path.basename(alignment)

    # Set TAPER output filename:
    expected_cleaned_alignment_file = (f'{output_folder}/'
                                       f'{re.sub(".fasta", ".cleaned.fasta", str(alignment_basename))}')

    command = (f'julia $(which correction_multi.jl) '
               f'-c {cleaning_cutoff} '
               f'-m - '
               f'-a N '
               f'{alignment}')

    logger.debug(f'TAPER command is: {command}')

    try:
        assert utils.file_exists_and_not_empty(expected_cleaned_alignment_file)
        logger.debug(f'Cleaned alignment exists for {alignment_basename}, skipping...')
        with lock:
            counter.value += 1

        return os.path.basename(expected_cleaned_alignment_file), None  # log list None as already processed

    except AssertionError:

        log_list_check_seqs = []

        try:
            result = subprocess.run(command, shell=True, universal_newlines=True, check=True, stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE)

            logger.debug(f'TAPER check_returncode() is: {result.check_returncode()}')
            logger.debug(f'TAPER stdout is: {result.stdout}')
            logger.debug(f'TAPER stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:

            logger.error(f'TAPER FAILED. Output is: {exc}')
            logger.error(f'TAPER stdout is: {exc.stdout}')
            logger.error(f'TAPER stderr is: {exc.stderr}')
            raise ValueError('There was an issue running TAPER. Check input files!')

        # Filter out empty sequences comprising only dashes, and post-TAPER alignments where all sequences
        # are either dashes or empty. If fewer than 4 'good' sequences are present, skip the gene:
        good_seqs = []
        seqs_all_dashes = []
        empty_seqs = []

        alignment_io = io.StringIO(result.stdout)
        seqs = SeqIO.parse(alignment_io, 'fasta')

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

            log_list_check_seqs.append(f'After running TAPER, the following sequences contained dashes, and have been '
                                       f'removed: {seqs_all_dashes_joined}')
        if empty_seqs:
            empty_seqs_joined = ', '.join(empty_seqs)

            log_list_check_seqs.append(f'After running TAPER, the following sequences were empty, and have been '
                                       f'removed: {empty_seqs_joined}')

        # Skip any filtered alignments with fewer than 4 sequences remaining:
        if len(good_seqs) < 4:

            log_list_check_seqs.append(f'{"[WARNING]:":10} After running TAPER, file {alignment_basename} contains '
                                       f'fewer than 4 good sequences, skipping gene!')
        else:
            with open(f'{expected_cleaned_alignment_file}', 'w') as filtered_taper_fasta:
                SeqIO.write(good_seqs, filtered_taper_fasta, 'fasta')

        with lock:
            counter.value += 1
            return os.path.basename(expected_cleaned_alignment_file), log_list_check_seqs

    finally:
        with lock:
            sys.stderr.write(f'\r{"[INFO]:":10} Finished generating cleaned alignment for file {alignment_basename}, '
                             f'{counter.value}/{num_files_to_process}')


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

    logger.info(f'\n{"[INFO]:":10} Generating alignments for fasta files using Clustal Omega...')
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

        for future in as_completed(future_results):

            try:
                check = future.result()

            except Exception as error:
                logger.error(f'Error raised: {error}')
                tb = traceback.format_exc()
                logger.error(f'traceback is:\n{tb}')
                sys.exit(1)

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
    Use Clustal Omega to align a fasta file of sequences, using the number of threads provided.
    Returns filename of the alignment produced.

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
        command = f'clustalo --in {fasta_file} --out {expected_alignment_file} --threads {threads} --verbose'
        logger.info(f'{"[INFO]:":10} Performing Clustal alignment with command: {command}')

        try:

            result = subprocess.run(command,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    check=True,
                                    shell=True)

            logger.debug(f'ClustalO check_returncode() is: {result.check_returncode()}')
            logger.debug(f'ClustalO stdout is: {result.stdout}')
            logger.debug(f'ClustalO stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'ClustalO FAILED. Output is: {exc}')
            logger.error(f'ClustalO stdout is: {exc.stdout}')
            logger.error(f'ClustalO stderr is: {exc.stderr}')
            raise ValueError('There was an issue running ClustalO. Check input files!')

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

    if not args.use_clustal:
        logger.debug(f'Running without no_stitched_contigs option - aligning with mafft only')

        alignments_output_folder = mafft_align_multiprocessing(
            gene_fasta_directory,
            algorithm=args.mafft_algorithm,
            adjust_direction=args.mafft_adjustdirection,
            pool_threads=args.pool,
            mafft_threads=args.threads,
            logger=logger)

        # Perform optional trimming with TrimAl:
        trimmed_output_folder = '03_alignments_trimmed'
        if not args.no_trimming:
            alignments_output_folder = run_trimal(alignments_output_folder,
                                                  trimmed_output_folder,
                                                  args_for_trimal_options=args,
                                                  logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping trimming step...')

        # Perform optional cleaning with TAPER (correction_multi.jl):
        if not args.no_cleaning:
            run_taper_multiprocessing(alignments_output_folder,
                                      no_trimming=args.no_trimming,
                                      cleaning_cutoff=args.cleaning_cutoff,
                                      pool_threads=args.pool,
                                      logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping cleaning step...')

    elif args.use_clustal:  # Align or re-align with Clustal Omega to better align short sequences
        logger.debug(f'Running with use_clustal option - aligning/realigning with Clustal Omega')

        if args.mafft_adjustdirection:  # Align with MAFFT first to revcomp any paralog sequences that need it:

            alignments_output_folder = mafft_align_multiprocessing(
                gene_fasta_directory,
                algorithm=args.mafft_algorithm,
                adjust_direction=args.mafft_adjustdirection,
                pool_threads=args.pool,
                mafft_threads=args.threads,
                logger=logger)

            alignments_output_folder = clustalo_align_multiprocessing(
                alignments_output_folder,
                pool_threads=args.pool,
                clustalo_threads=args.threads,
                logger=logger)

        else:
            alignments_output_folder = clustalo_align_multiprocessing(
                gene_fasta_directory,
                pool_threads=args.pool,
                clustalo_threads=args.threads,
                logger=logger)

        # Perform optional trimming with TrimAl:
        trimmed_output_folder = '03_alignments_trimmed'
        if not args.no_trimming:
            alignments_output_folder = run_trimal(alignments_output_folder,
                                                  trimmed_output_folder,
                                                  args_for_trimal_options=args,
                                                  logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping trimming step...')

        # Perform optional cleaning with HmmCleaner.pl
        if not args.no_cleaning:
            run_taper_multiprocessing(alignments_output_folder,
                                      no_trimming=args.no_trimming,
                                      pool_threads=args.pool,
                                      logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping cleaning step...')

    logger.info(f'\n{"[INFO]:":10} Finished aligning and trimming/cleaning input files.')

