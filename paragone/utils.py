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
import io
import pstats
import cProfile
import re
from textwrap import TextWrapper


def check_inputs(directory_suffix_dict,
                 file_list,
                 logger=None):
    """
    Checks that provided directories and files both exist and are not empty.

    :param directory_suffix_dict:
    :param file_list:
    :param logging.Logger logger: a logger object
    :return:
    """

    # Check that inpout directories exist, contains files with the expected suffix, and the files are not empty:
    for directory, expected_file_suffix in directory_suffix_dict.items():
        logger.debug(f'Checking directory {directory} for files with suffix {expected_file_suffix}')
        if not os.path.isdir(directory):
            logger.error(f'{"[ERROR]:":10} Directory not found: {directory}')
            sys.exit(1)
        else:
            logger.debug(f'Directory {directory} exists, proceeding...')

        expected_files = glob.glob(f'{directory}/*{expected_file_suffix}')
        if not expected_files:
            logger.error(f'{"[ERROR]:":10} Directory "{directory}" contains no files with suffix:'
                         f' {expected_file_suffix}')
            sys.exit(1)
        else:
            logger.debug(f'Expected_files are {expected_files}, proceeding...')

        empty_files = []
        for item in expected_files:
            if not file_exists_and_not_empty(item):
                empty_files.append(os.path.basename(item))
            else:
                logger.debug(f'Expected file {item} is not empty, proceeding...')

        if empty_files:
            joined = ', '.join(empty_files)
            logger.error(f'{"[ERROR]:":10} The following file(s) in directory "{directory}" are empty: {joined}')
            sys.exit(1)

    # Check that any files provided exist, and are not empty:
    missing_files = []
    empty_files = []
    for item in file_list:
        if not os.path.isfile(item):
            missing_files.append(item)
        else:
            logger.debug(f'Expected file {item} exists, proceeding...')

    for item in file_list:
        if item not in missing_files and not file_exists_and_not_empty(item):
            empty_files.append(item)
        else:
            logger.debug(f'Expected file {item} is not empty, proceeding...')

    if missing_files:
        joined = ', '.join(missing_files)
        logger.error(f'{"[ERROR]:":10} The following file(s) do not exist: {joined}')
        sys.exit(1)

    if empty_files:
        joined = ', '.join(empty_files)
        logger.error(f'{"[ERROR]:":10} The following file(s) are empty: {joined}')
        sys.exit(1)


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

    logger.info(f'\n{"[INFO]:":10} Resolving any polytomies in trees...')

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


def parse_ingroup_and_outgroup_file(in_out_file, logger=None):
    """
    Parse an input text file and return a list of ingroup taxa and a list of outgroup taxa.

    :param str in_out_file: path to the text file containing ingroup and outgroup designations
    :param logging.Logger logger: a logger object
    :return list ingroups, outgroups: lists of ingroup taxa and outgroup taxa
    """

    ingroups = []
    outgroups = []

    with open(in_out_file, 'r') as in_out_handle:
        for line in in_out_handle:
            if len(line) < 3:
                logger.debug(f'Skipping line {line} in in_out_file {os.path.basename(in_out_file)} as len < 3')
                continue
            sample = line.strip().split("\t")
            if sample[0] == "IN":
                ingroups.append(sample[1])
            elif sample[0] == "OUT":
                outgroups.append(sample[1])
            else:
                logger.error(f'{"[ERROR]:":10} Check in_and_outgroup_list file format for the following line:')
                logger.error(f'\n{" " * 10} {line}')
                sys.exit(1)

    # Check if there are taxa designated as both ingroup AND outgroup:
    if len(set(ingroups) & set(outgroups)) > 0:
        logger.error(f'{"[ERROR]:":10} Taxon ID {set(ingroups)} & {set(outgroups)} are in both ingroup and outgroup!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} There are {len(ingroups)} ingroup taxa and {len(outgroups)} outgroup taxa on the'
                f' {os.path.basename(in_out_file)} file provided')

    return ingroups, outgroups


def cprofile_to_csv(profile_binary_file):
    """
    Takes a cProfile.Profile object, converts it to human-readable format, and returns the data in *.csv format for
    writing to file.

    From: https://gist.github.com/ralfstx/a173a7e4c37afa105a66f371a09aa83e

    :param cProfile.Profile profile_binary_file: a cProfile.Profile object
    :return str: human-readable data from the cProfile run, in *.csv format
    """

    out_stream = io.StringIO()
    pstats.Stats(profile_binary_file, stream=out_stream).sort_stats('cumtime').print_stats()
    result = out_stream.getvalue()
    result = 'ncalls' + result.split('ncalls')[-1]  # chop off header lines
    lines = [','.join(line.rstrip().split(None, 5)) for line in result.split('\n')]

    return '\n'.join(lines)


def fill_forward_slash(text, width=70, **kwargs):
    """
    Fill a single paragraph of text, returning a new string.

    Reformat the single paragraph in 'text' to fit in lines of no more
    than 'width' columns, and return a new string containing the entire
    wrapped paragraph.  As with wrap(), tabs are expanded and other
    whitespace characters converted to space.  See TextWrapper class for
    available keyword args to customize wrapping behaviour.

    This function uses the subclass TextWrapperForwardSlash.
    """
    w = TextWrapperForwardSlash(width=width, **kwargs)
    return w.fill(text)


class TextWrapperForwardSlash(TextWrapper):
    """
    Subclasses textwrap.TextWrapper and alters regex so that break_on_hyphens corresponds to forward slashes rather
    than hyphens. Used for wrapping long path strings.

    Change: letter = r'[^\d\W]' -> letter = r'[\w-]'
    """

    _whitespace = '\t\n\x0b\x0c\r '
    word_punct = r'[\w!"\'&.,?]'
    letter = r'[\w-]'
    whitespace = r'[%s]' % re.escape(_whitespace)
    nowhitespace = '[^' + whitespace[1:]

    wordsep_re = re.compile(r'''
        ( # any whitespace
          %(ws)s+
        | # em-dash between words
          (?<=%(wp)s) -{2,} (?=\w)
        | # word, possibly hyphenated
          %(nws)s+? (?:
            # hyphenated word
              (?: (?<=%(lt)s{2}/) | (?<=%(lt)s/%(lt)s/))
              (?= %(lt)s /? %(lt)s)
            | # end of word
              (?=%(ws)s|\Z)
            | # em-dash
              (?<=%(wp)s) (?=-{2,}\w)
            )
        )''' % {'wp': word_punct,
                'lt': letter,
                'ws': whitespace,
                'nws': nowhitespace},
        re.VERBOSE)
    del word_punct, letter, nowhitespace

    def __init__(self,
                 width=70,
                 initial_indent="",
                 subsequent_indent="",
                 expand_tabs=True,
                 replace_whitespace=True,
                 fix_sentence_endings=False,
                 break_long_words=True,
                 drop_whitespace=True,
                 break_on_forward_slash=True,
                 tabsize=8,
                 *,
                 max_lines=None,
                 placeholder=' [...]'):
        super().__init__(width=width,
                         initial_indent=initial_indent,
                         subsequent_indent=subsequent_indent,
                         expand_tabs=expand_tabs,
                         replace_whitespace=replace_whitespace,
                         fix_sentence_endings=fix_sentence_endings,
                         break_long_words=break_long_words,
                         drop_whitespace=drop_whitespace,
                         break_on_hyphens=break_on_forward_slash,
                         tabsize=tabsize,
                         max_lines=max_lines,
                         placeholder=placeholder)