#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au

"""
This module contains some general functions and classes.
"""

import os
import sys
import shutil
import logging
import datetime
import glob
from ete3 import Tree


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


def createfolder(directory):
    """
    Attempts to create a directory named after the name provided, and provides an error message on failure
    """

    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print(f'{"[ERROR]:":10} Error creating directory: {directory}')
        sys.exit(1)


def file_exists_and_not_empty(file_name):
    """
    Check if file exists and is not empty by confirming that its size is not 0 bytes
    """

    return os.path.isfile(file_name) and not os.path.getsize(file_name) == 0


def done_callback(future_returned):
    """
    Callback function for ProcessPoolExecutor futures; gets called when a future is cancelled or 'done'.
    """
    if future_returned.cancelled():
        print(f'{future_returned}: cancelled')
        return
    elif future_returned.done():
        error = future_returned.exception()
        if error:
            print(f'{future_returned}: error returned: {error}')
        else:
            result = future_returned.result()
            return result


def check_dependencies(logger=None):
    """
    Checks for the presence of executables. Returns a boolean.

    :param None, logging.Logger logger: a logger object
    return: bool: True if all dependencies are found, else False
    """

    executables = ['mafft',
                   'clustalo',
                   'iqtree',
                   'FastTreeMP',
                   'HmmCleaner.pl']

    logger.info(f'{"[INFO]:":10} Checking for external dependencies:\n')

    all_executable_found = True
    for executable in executables:
        executable_loc = shutil.which(executable)
        if executable:
            logger.info(f'{executable:20} found at {executable_loc}')
        else:
            logger.info(f'{executable:20} not found in your $PATH!')
            all_executable_found = False

    logger.info('')

    return all_executable_found


def resolve_polytomies(treefile_directory, logger=None):
    """
    Iterates over tree newick files in a directory. For each tree, checks for polyotomies and arbitrarily resolves
    them. Overwrites the original tree with the resolved tree.

    :param str treefile_directory: name of directory containing tree files in newick format
    :param logging.Logger logger: a logger object
    :return:
    """

    for treefile in glob.glob(f'{treefile_directory}/*.treefile'):
        tree_basename = os.path.basename(treefile)
        logger.debug(f'Examining tree {tree_basename} for polytomies...')

        # Randomly resolve any polytomies (as FastTree can generate polytomies):
        ete3_tree = Tree(treefile)
        polytomy = False

        for node in ete3_tree.iter_descendants():
            if len(node.children) > 2:
                logger.debug(f'Tree {tree_basename} has a polytomy at node {node.write()}')
                polytomy = True

        if polytomy:
            logger.debug(f'Polytomies in tree {tree_basename} will be randomly resolved! The original tree is: '
                         f'{ete3_tree.write()}')

            ete3_tree.resolve_polytomy(recursive=True, default_dist='0.01')
            for node in ete3_tree.iter_descendants():
                if len(node.children) > 2:
                    raise ValueError(f'Tree {tree_basename} still has a polytomy at node {node.write()}!')
            ete3_tree.unroot()
            ete3_tree.write(outfile=f'{treefile}', format=0)

        else:
            logger.debug(f'Tree {tree_basename} has no polytomies - keeping original file.')

