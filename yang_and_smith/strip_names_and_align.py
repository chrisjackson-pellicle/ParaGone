#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
Strip sequence names, ready for concatenation, and aligns each fasta file using mafft/clustalo.
"""

import sys
import os
import re
import subprocess
import textwrap
import glob
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline, ClustalOmegaCommandline, MuscleCommandline
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait

from yang_and_smith import utils


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
                mafft_cline = (MafftCommandline(auto='true', thread=threads, input=fasta_file))
            else:
                mafft_cline = (MafftCommandline(algorithm, thread=threads, input=fasta_file))

            logger.info(f'{"[INFO]:":10} Performing MAFFT alignment with command: {mafft_cline}')
            stdout, stderr = mafft_cline()

            logger.debug(f'stdout is: {stdout}')
            logger.debug(f'stderr is: {stderr}')

            with open(expected_alignment_file, 'w') as alignment_file:
                alignment_file.write(stdout)

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
                raise ValueError('There was an issue running trimal. Check input files!')

        with lock:
            counter.value += 1
            logger.debug(f'Aligned file {fasta_file_basename}')

        return os.path.basename(expected_alignment_file)

    finally:
        with lock:
            sys.stderr.write(f'\r{"[INFO]:":10} Finished generating alignment for file {fasta_file_basename},'
                             f' {counter.value}/{num_files_to_process}')


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
    output_folder = f'24_{input_folder_basename.lstrip("_23")}_alignments'
    utils.createfolder(output_folder)

    if use_muscle:
        logger.info(f'{"[INFO]:":10} Generating alignments for fasta files using MUSCLE...')
    else:
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
        future_results = [pool.submit(mafft_or_muscle_align,
                                      fasta_file,
                                      algorithm,
                                      output_folder,
                                      counter,
                                      lock,
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
            raise ValueError('There was an issue running trimal. Check input files!')

        with lock:
            counter.value += 1
            logger.debug(f'Aligned file {fasta_file_basename}')

        return os.path.basename(expected_alignment_file)

    finally:
        with lock:
            sys.stderr.write(f'\r{"[INFO]:":10} Finished generating alignment for file {fasta_file_basename},'
                             f' {counter.value}/{num_files_to_process}')


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


def strip_names_for_concat(selected_alignment_directory):
    """
    Strips everything after a dot ('.') from the name of each sequence in an alignment file.
    Returns the name of the output folder.

    :param str selected_alignment_directory: path to input folder containing fasta alignment files
    :return str: path to the output folder produced
    """

    input_folder_basename = os.path.basename(selected_alignment_directory)
    output_folder = f'23_{input_folder_basename}_stripped_names'
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

    return output_folder


def main(args):
    """
    Entry point for the resolve_paralogs.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = utils.setup_logger(__name__, '00_logs_and_reports_resolve_paralogs/logs/15_strip_names_and_align')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand strip_names_and_align was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')
    logger.debug(args)

    # Checking input directories and files:
    directory_suffix_dict = {args.selected_alignment_directory: '.fasta'}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Create output folder for pruned trees:
    # output_folder = f'{os.path.basename(args.selected_alignment_directory)}_stripped_names'
    # utils.createfolder(output_folder)

    stripped_names_output_folder = strip_names_for_concat(args.selected_alignment_directory)

    if not args.no_stitched_contigs:  # i.e. if it's a standard run with stitched contigs produced.
        logger.debug(f'Running without no_stitched_contigs option - aligning with mafft or muscle only')
        alignments_output_folder = mafft_or_muscle_align_multiprocessing(
            stripped_names_output_folder,
            algorithm=args.mafft_algorithm,
            pool_threads=args.pool,
            mafft_threads=args.threads,
            no_stitched_contigs=args.no_stitched_contigs,
            use_muscle=args.use_muscle,
            logger=logger)

    elif args.no_stitched_contigs:  # Re-align with Clustal Omega.
        logger.debug(f'Running with no_stitched_contigs option - realigning with clustal omega')
        alignments_output_folder = mafft_or_muscle_align_multiprocessing(
            stripped_names_output_folder,
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

    logger.info(f'\n{"[INFO]:":10} Finished stripping fasta sequence names and final alignments.')
