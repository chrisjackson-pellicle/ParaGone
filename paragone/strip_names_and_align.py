#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
Strip sequence names, ready for concatenation, and aligns each fasta file using mafft/clustalo.
"""

import sys
import os
import re
import textwrap
import glob
from Bio import SeqIO, AlignIO
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait, as_completed
import traceback
import subprocess

from paragone import utils
from paragone.align_and_clean import run_trimal


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

        return os.path.basename(expected_alignment_file)

    except AssertionError:

        if algorithm == 'auto':
            command = f'mafft --auto --thread {threads} {fasta_file} > {expected_alignment_file}'
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

        with lock:
            counter.value += 1
            logger.debug(f'Aligned file {fasta_file_basename}')

        return os.path.basename(expected_alignment_file)

    finally:
        with lock:
            sys.stderr.write(f'\r{"[INFO]:":10} Finished generating alignment for file {fasta_file_basename},'
                             f' {counter.value}/{num_files_to_process}')


def mafft_align_multiprocessing(fasta_to_align_folder,
                                alignments_output_folder,
                                algorithm='auto',
                                pool_threads=1,
                                mafft_threads=2,
                                logger=None):
    """
    Generate alignments via function <mafft_align> using multiprocessing.

    :param str fasta_to_align_folder: path to folder containing input fasta files
    :param str alignments_output_folder: name of folder for output alignments
    :param str algorithm: algorithm to use for mafft alignment; default is 'auto'
    :param int pool_threads: number of alignments to run concurrently
    :param int mafft_threads: number of threads to use for each concurrent alignment
    :param logging.Logger logger: a logger object
    :return str output_folder: name of the output folder containing alignments
    """

    utils.createfolder(alignments_output_folder)

    logger.info(f'{"[INFO]:":10} Generating alignments for fasta files using MAFFT...')

    # Filter out any input files with fewer than four sequences:
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
        future_results = [pool.submit(mafft_align,
                                      fasta_file,
                                      algorithm,
                                      alignments_output_folder,
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

    alignment_list = [alignment for alignment in glob.glob(f'{alignments_output_folder}/*.aln.fasta') if
                      utils.file_exists_and_not_empty(alignment)]

    logger.debug(f'{len(alignment_list)} alignments generated from {len(future_results)} fasta files...')

    return alignments_output_folder


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


def clustalo_align_multiprocessing(fasta_to_align_folder,
                                   alignments_output_folder,
                                   pool_threads=1,
                                   clustalo_threads=1,
                                   logger=None):
    """
    Generate alignments via function <clustalo_align> using multiprocessing.

    :param fasta_to_align_folder: path to folder containing input fasta alignment files
    :param str alignments_output_folder: name of folder for output alignments
    :param int pool_threads: number of alignments to run concurrently
    :param int clustalo_threads: number of threads to use for each concurrent alignment
    :param logging.Logger logger: a logger object
    :return str output_folder: name of the output folder containing alignments
    """

    utils.createfolder(alignments_output_folder)

    logger.info(f'{"[INFO]:":10} Generating alignments for fasta files using Clustal Omega...')

    target_genes = [file for file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(clustalo_align,
                                      fasta_file,
                                      alignments_output_folder,
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

    alignment_list = [alignment for alignment in glob.glob(f'{alignments_output_folder}/*.aln.fasta') if
                      utils.file_exists_and_not_empty(alignment)]

    logger.debug(f'{len(alignment_list)} alignments generated from {len(future_results)} fasta files...')

    return alignments_output_folder


def strip_names_for_concat(selected_alignment_directory,
                           output_folder):
    """
    Strips everything after a dot ('.') from the name of each sequence in an alignment file.
    Returns the name of the output folder.

    :param str selected_alignment_directory: path to input folder containing fasta alignment files
    :param str output_folder: path to output folder for sequences with stripped names
    # :return str: path to the output folder produced
    :return
    """

    utils.createfolder(output_folder)

    for alignment in glob.glob(f'{selected_alignment_directory}/*.fasta'):
        alignment_basename = os.path.basename(alignment)
        seqs = AlignIO.read(alignment, 'fasta')
        for seq in seqs:
            seq.name = f'{seq.name.split(".")[0]}'
            seq.id = f'{seq.id.split(".")[0]}'
            seq.description = ''
        with open(f'{output_folder}/{re.sub(".fasta", "_stripped.fasta", alignment_basename)}', 'w') as stripped:
            AlignIO.write(seqs, stripped, 'fasta')


def main(args,
         report_directory,
         selected_alignment_directory=None,
         logger=None):
    """
    Entry point for the paragone_main.py script

    :param args: argparse namespace with subparser options for function main()
    :param str report_directory: path to directory for report files
    :param str selected_alignment_directory: name of directory with selected sequences
    :param logging.Logger logger: a logger object
    :return:
    """

    logger.debug(f'{"[INFO]:":10} Module strip_names_and_align was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]),
                         width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11, break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    logger.info('')
    logger.info(f'{"[INFO]:":10} ======> STRIP PARALOG NAME SUFFIX AND ALIGN <======\n')

    # Checking input directories and files:
    directory_suffix_dict = {selected_alignment_directory: '.fasta'}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    if re.search('MO', selected_alignment_directory):
        stripped_names_folder = f'20_MO_stripped_names'
        alignments_output_folder = f'23_MO_final_alignments'
        alignments_output_folder_trimmed = f'26_MO_final_alignments_trimmed'
    if re.search('MI', selected_alignment_directory):
        stripped_names_folder = f'21_MI_stripped_names'
        alignments_output_folder = f'24_MI_final_alignments'
        alignments_output_folder_trimmed = f'27_MI_final_alignments_trimmed'
    if re.search('RT', selected_alignment_directory):
        stripped_names_folder = f'22_RT_stripped_names'
        alignments_output_folder = f'25_RT_final_alignments'
        alignments_output_folder_trimmed = f'28_RT_final_alignments_trimmed'

    strip_names_for_concat(selected_alignment_directory,
                           stripped_names_folder)

    if not args.use_clustal:
        logger.debug(f'Running without --use_clustal option - aligning with MAFFT only')

        alignments_output_folder = mafft_align_multiprocessing(
            stripped_names_folder,
            alignments_output_folder,
            algorithm=args.mafft_algorithm,
            pool_threads=args.pool,
            mafft_threads=args.threads,
            logger=logger)

        # Perform optional trimming with TrimAl:
        if not args.no_trimming:
            alignments_output_folder = run_trimal(alignments_output_folder,
                                                  alignments_output_folder_trimmed,
                                                  args_for_trimal_options=args,
                                                  logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping trimming step...')

    elif args.use_clustal:  # Align with Clustal Omega.
        logger.debug(f'Running with --use_clustal option - aligning with Clustal Omega only')

        alignments_output_folder = clustalo_align_multiprocessing(
            stripped_names_folder,
            alignments_output_folder,
            pool_threads=args.pool,
            clustalo_threads=args.threads,
            logger=logger)

        # Perform optional trimming with TrimAl:
        if not args.no_trimming:
            alignments_output_folder = run_trimal(alignments_output_folder,
                                                  alignments_output_folder_trimmed,
                                                  args_for_trimal_options=args,
                                                  logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping trimming step...')

    logger.info(f'\n{"[INFO]:":10} Finished stripping fasta sequence names and final alignments.')
