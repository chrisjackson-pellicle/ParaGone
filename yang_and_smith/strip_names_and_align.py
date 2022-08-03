#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au

"""
Strip names ready for concatentation and aligns each fasta file using mafft/clustalo.

"""

import sys
import os
import re
import subprocess
import shutil
from collections import defaultdict
import textwrap
import glob
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline, ClustalOmegaCommandline, MuscleCommandline
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait

from yang_and_smith import utils


def mafft_align(fasta_file, algorithm, output_folder, counter, lock, num_files_to_process, threads=2,
                no_supercontigs=False, use_muscle=False):
    """
    Uses mafft to align a fasta file of sequences, using the algorithm and number of threads provided. Returns filename
    of the alignment produced.
    """

    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'

    try:
        assert file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
        with lock:
            counter.value += 1
        return os.path.basename(expected_alignment_file)
    except AssertionError:
        if use_muscle:
            # print(fasta_file_basename)
            logger.info('Alignment will be performed using MUSCLE rather than MAFFT!')
            muscle_cline = MuscleCommandline(input=fasta_file, out=expected_alignment_file)
            stdout, stderr = muscle_cline()
        else:
            if algorithm == 'auto':
                mafft_cline = (MafftCommandline(auto='true', thread=threads, input=fasta_file))
            else:
                mafft_cline = (MafftCommandline(algorithm, thread=threads, input=fasta_file))
            logger.info(f'Performing MAFFT alignment with command" {mafft_cline}')
            stdout, stderr = mafft_cline()
            with open(expected_alignment_file, 'w') as alignment_file:
                alignment_file.write(stdout)
        with lock:
            counter.value += 1
        logger.debug(f'Aligned file {fasta_file_basename}')
        return os.path.basename(expected_alignment_file)
    finally:
        print(f'\rFinished generating alignment for file {fasta_file_basename}, '
              f'{counter.value}/{num_files_to_process}', end='')


def mafft_align_multiprocessing(fasta_to_align_folder, algorithm='linsi', pool_threads=1,
                                mafft_threads=2, no_supercontigs=False, use_muscle=False):
    """
    Generate alignments via function <align_targets> using multiprocessing.
    """

    input_folder_basename = os.path.basename(fasta_to_align_folder)
    output_folder = f'{input_folder_basename}_alignments'
    # print(f'output_folder from mafft_align_multiprocessing: {output_folder}')
    createfolder(output_folder)

    logger.info(f'Generating alignments for fasta files in folder {fasta_to_align_folder}...\n')
    target_genes = [file for file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta'))]

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
                                      no_supercontigs=no_supercontigs,
                                      use_muscle=use_muscle)
                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                      file_exists_and_not_empty(alignment)]
    logger.info(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')

    return output_folder


def clustalo_align(fasta_file, output_folder, counter, lock, num_files_to_process, threads=2):
    """
    Uses clustal omega to align a fasta file of sequences, using the algorithm and number of threads provided. Returns
    filename
    of the alignment produced.
    """

    createfolder(output_folder)
    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{fasta_file_basename}'

    try:
        assert file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
    except AssertionError:
        # clustalomega_cline = ClustalOmegaCommandline(infile=fasta_file, outfile=expected_alignment_file,
        #                                              verbose=True, auto=True, threads=threads, pileup=True)
        clustalomega_cline = ClustalOmegaCommandline(infile=fasta_file, outfile=expected_alignment_file,
                                                     verbose=True, auto=True, threads=threads)
        clustalomega_cline()
    finally:
        with lock:
            counter.value += 1
            print(f'\rFinished generating alignment for file {fasta_file_basename}, '
                  f'{counter.value}/{num_files_to_process}', end='')
        return os.path.basename(expected_alignment_file)


def clustalo_align_multiprocessing(fasta_to_align_folder, pool_threads=1, clustalo_threads=2):
    """
    Generate alignments via function <clustalo_align> using multiprocessing.
    """

    input_folder_basename = os.path.basename(fasta_to_align_folder)
    output_folder = f'{input_folder_basename}_clustal'
    # print(f'output_folder from clustal_align_multiprocessing: {output_folder}')
    createfolder(output_folder)

    logger.info('Generating alignments for fasta files using clustal omega...\n')
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
                                      threads=clustalo_threads)
                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                      file_exists_and_not_empty(alignment)]
    logger.info(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')

    return output_folder


def strip_names_for_concat(input_folder):
    """
    Strips everything after a dot ('.') from the name of each sequence in an alignment file.
    Returns the name of the output folder.
    """

    input_folder_basename = os.path.basename(input_folder)
    output_folder = f'{input_folder_basename}_stripped_names'
    # print(f'output_folder from mafft_align_multiprocessing: {output_folder}')
    createfolder(output_folder)

    for alignment in glob.glob(f'{input_folder}/*.fa'):
        alignment_basename = os.path.basename(alignment)
        seqs = AlignIO.read(alignment, 'fasta')
        for seq in seqs:
            seq.name = f'{seq.name.split(".")[0]}'
            seq.id = f'{seq.id.split(".")[0]}'
            seq.description = ''
        with open(f'{output_folder}/{re.sub(".fa", "_stripped.fasta", alignment_basename)}', 'w') as stripped:
            AlignIO.write(seqs, stripped, 'fasta')
    return output_folder


def main(args):
    """
    Entry point for the resolve_paralogs.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = utils.setup_logger(__name__, 'logs_resolve_paralogs/13_strip_names_and_align')

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

    # Create output folder for pruned trees:
    output_folder = f'{os.path.basename(args.selected_alignment_directory)}_stripped_names'
    utils.createfolder(output_folder)

    # output_folder = strip_names_for_concat(args.selected_alignment_directory)

    # if not results.no_supercontigs:  # i.e. if it's a standard run.
    #     mafft_align_multiprocessing(output_folder,
    #                                 algorithm=results.mafft_algorithm,
    #                                 pool_threads=results.threads_pool,
    #                                 mafft_threads=results.threads_mafft,
    #                                 no_supercontigs=results.no_supercontigs,
    #                                 use_muscle=results.use_muscle)
    #
    # elif results.no_supercontigs:  # re-align with Clustal Omega.
    #     alignments_output_folder = mafft_align_multiprocessing(output_folder,
    #                                                            algorithm=results.mafft_algorithm,
    #                                                            pool_threads=results.threads_pool,
    #                                                            mafft_threads=results.threads_mafft,
    #                                                            no_supercontigs=results.no_supercontigs,
    #                                                            use_muscle=results.use_muscle)
    #
    #     clustalo_align_multiprocessing(alignments_output_folder,
    #                                    pool_threads=results.threads_pool,
    #                                    clustalo_threads=results.threads_mafft)


