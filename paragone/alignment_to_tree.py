#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
Takes a trimmed, Hmmcleaned alignment, and produces a tree via FastTreeMP or IQTree
"""

import logging
import sys
import os
import glob
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait
import subprocess
import textwrap

from paragone import utils


def fasttree_multiprocessing(alignments_folder,
                             pool=1,
                             threads=1,
                             bootstraps=False,
                             logger=None):
    """
    Generate FastTree trees using multiprocessing.

    :param str alignments_folder: name of folder containing alignments
    :param int pool: number of concurrent trees to run; default is 1
    :param int threads: number threads for each concurrent tree; default is 1
    :param bool bootstraps: if False, turn off support calculation
    :param logging.Logger logger: a logger object
    :return str output_folder: path to output folder containing tree files
    """

    # input_folder_basename = os.path.basename(alignments_folder)
    output_folder = f'05_trees_pre_quality_control'
    utils.createfolder(output_folder)

    fill = textwrap.fill(f'{"[INFO]:":10} Generating trees from alignments using FastTreeMP. Tree files will be '
                         f'written to directory: "{output_folder}".',
                         width=90, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)

    logger.info(f'{fill}')

    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(fasttree,
                                      alignment,
                                      output_folder,
                                      threads,
                                      counter,
                                      lock,
                                      num_files_to_process=len(alignments),
                                      bootstraps=bootstraps,
                                      logger=logger)
                          for alignment in alignments]
        for future in future_results:
            future.add_done_callback(utils.done_callback)

        wait(future_results, return_when="ALL_COMPLETED")

    tree_list = [tree for tree in glob.glob(f'{output_folder}/*.treefile') if utils.file_exists_and_not_empty(tree)]

    logger.debug(f'{len(tree_list)} trees generated from {len(future_results)} fasta files...')

    return output_folder


def fasttree(alignment_file,
             output_folder,
             threads,
             counter,
             lock,
             num_files_to_process,
             bootstraps=False,
             logger=None):
    """
    Generate trees from alignments using FastTreeMP

    :param str alignment_file: path to alignment file
    :param str output_folder: name of output folder for tree
    :param int threads: number of threads to use for FastTreeMP
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int num_files_to_process: total number of alignment fasta files for tree generation
    :param bool bootstraps: if False, turn off support calculation
    :param logging.Logger logger: a logger object
    :return str expected_output_file: name of expected output tree file
    """

    alignment_file_basename = os.path.basename(alignment_file)
    expected_output_file = f'{output_folder}/{alignment_file_basename}.treefile'

    try:
        assert utils.file_exists_and_not_empty(expected_output_file)
        logger.debug(f'Output exists for {expected_output_file}, skipping...')
        with lock:
            counter.value += 1

        return os.path.basename(expected_output_file)

    except AssertionError:
        try:
            if bootstraps:
                fasttree_command = f'export OMP_NUM_THREADS={threads}; ' \
                                   f'FastTreeMP -gtr -nt < {alignment_file} > {expected_output_file}'
                result = subprocess.run(fasttree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
                logger.debug(f'FastTreeMP check_returncode() is: {result.check_returncode()}')
                logger.debug(f'FastTreeMP stdout is: {result.stdout}')
                logger.debug(f'FastTreeMP stderr is: {result.stderr}')

            else:
                fasttree_command = f'export OMP_NUM_THREADS={threads};' \
                                   f'FastTreeMP -gtr -nt -nosupport < {alignment_file} > {expected_output_file}'
                result = subprocess.run(fasttree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
                logger.debug(f'FastTreeMP check_returncode() is: {result.check_returncode()}')
                logger.debug(f'FastTreeMP stdout is: {result.stdout}')
                logger.debug(f'FastTreeMP stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'FastTreeMP FAILED. Output is: {exc}')
            logger.error(f'FastTreeMP stdout is: {exc.stdout}')
            logger.error(f'FastTreeMP stderr is: {exc.stderr}')

        with lock:
            counter.value += 1

        return os.path.basename(expected_output_file)

    finally:
        sys.stderr.write(f'\r{"[INFO]:":10} Finished generating tree {os.path.basename(expected_output_file)},'
                         f' {counter.value}/{num_files_to_process}')


def iqtree_multiprocessing(alignments_folder,
                           pool=1,
                           threads=2,
                           bootstraps=False,
                           logger=None):
    """
    Generate IQTree trees using multiprocessing.

    :param str alignments_folder: name of folder containing alignments
    :param int pool: number of concurrent trees to run; default is 1
    :param int threads: number threads for each concurrent tree; default is 1
    :param bool bootstraps: if False, turn off support calculation
    :param logging.Logger logger: a logger object
    :return str output_folder: path to output folder containing tree files
    """

    # input_folder_basename = os.path.basename(alignments_folder)
    output_folder = f'05_trees_pre_quality_control'
    utils.createfolder(output_folder)

    fill = textwrap.fill(f'{"[INFO]:":10} Generating trees from alignments using IQTREE. Tree files will be '
                         f'written to directory: "{output_folder}".',
                         width=90, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)

    logger.info(f'{fill}')

    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(iqtree,
                                      alignment,
                                      output_folder,
                                      threads,
                                      counter,
                                      lock,
                                      num_files_to_process=len(alignments),
                                      bootstraps=bootstraps,
                                      logger=logger)
                          for alignment in alignments]
        for future in future_results:
            future.add_done_callback(utils.done_callback)

        wait(future_results, return_when="ALL_COMPLETED")

    tree_list = [tree for tree in glob.glob(f'{output_folder}/*.treefile') if utils.file_exists_and_not_empty(tree)]

    logger.debug(f'{len(tree_list)} trees generated from {len(future_results)} fasta files...')

    return output_folder


def iqtree(alignment_file,
           output_folder,
           threads,
           counter,
           lock,
           num_files_to_process,
           bootstraps=False,
           logger=None):
    """
    Generate trees from alignments using IQTREE

    :param str alignment_file: path to alignment file
    :param str output_folder: name of output folder for tree
    :param int threads: number of threads to use for FastTreeMP
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int num_files_to_process: total number of alignment fasta files for tree generation
    :param bool bootstraps: if False, turn off support calculation
    :param logging.Logger logger: a logger object
    :return str expected_output_file: name of expected output tree file
    """

    alignment_file_basename = os.path.basename(alignment_file)
    expected_output_file = f'{output_folder}/{alignment_file_basename}.treefile'

    try:
        assert utils.file_exists_and_not_empty(expected_output_file)
        logger.debug(f'Output exists for {expected_output_file}, skipping...')
        with lock:
            counter.value += 1

        return os.path.basename(expected_output_file)

    except AssertionError:
        try:
            if bootstraps:
                iqtree_command = f'iqtree -redo -pre {output_folder}/{alignment_file_basename} -s {alignment_file} ' \
                                 f'-m GTR+G -bb 1000 -bnni -nt {str(threads)} -quiet'
                result = subprocess.run(iqtree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
                logger.debug(f'IQTREE check_returncode() is: {result.check_returncode()}')
                logger.debug(f'IQTREE stdout is: {result.stdout}')
                logger.debug(f'IQTREE stderr is: {result.stderr}')

            else:
                iqtree_command = f'iqtree -redo -pre {output_folder}/{alignment_file_basename} -s {alignment_file} ' \
                                 f'-m GTR+G -nt {str(threads)} -quiet'
                result = subprocess.run(iqtree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
                logger.debug(f'IQTREE check_returncode() is: {result.check_returncode()}')
                logger.debug(f'IQTREE stdout is: {result.stdout}')
                logger.debug(f'IQTREE stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'IQTREE FAILED. Output is: {exc}')
            logger.error(f'IQTREE stdout is: {exc.stdout}')
            logger.error(f'IQTREE stderr is: {exc.stderr}')

        with lock:
            counter.value += 1

        return os.path.basename(expected_output_file)

    finally:
        sys.stderr.write(f'\r{"[INFO]:":10} Finished generating output {os.path.basename(expected_output_file)},'
                         f' {counter.value}/{num_files_to_process}')


def main(args, logger=None):
    """
    Entry point for the paragone_main.py script

    :param args: argparse namespace with subparser options for function main()
    :param logging.Logger logger: a logger object
    :return:
    """

    logger.debug(f'{"[INFO]:":10} Module alignment_to_tree was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    logger.info(f'{"[INFO]:":10} ======> GENERATING TREES FROM PARALOG ALIGNMENTS <======\n')

    # Checking input directories and files:
    directory_suffix_dict = {args.alignment_directory: '.fasta'}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    if args.use_fasttree:
        trees_folder = fasttree_multiprocessing(args.alignment_directory,
                                                pool=args.pool,
                                                threads=args.threads,
                                                bootstraps=args.generate_bootstraps,
                                                logger=logger)

        utils.resolve_polytomies(trees_folder, logger=logger)

    else:
        trees_folder = iqtree_multiprocessing(args.alignment_directory,
                                              pool=args.pool,
                                              threads=args.threads,
                                              bootstraps=args.generate_bootstraps,
                                              logger=logger)

        utils.resolve_polytomies(trees_folder, logger=logger)

    logger.info(f'{"[INFO]:":10} Finished generating trees from alignments.')

