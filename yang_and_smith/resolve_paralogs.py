#!/usr/bin/env python

"""
Adapted Yang and Smith paralogy resolution pipeline version 1.0.0 release candidate (July 2022)

Author: Chris Jackson chris.jackson@rbg.vic.gov.au

"""

import argparse
import os
import sys
import shutil
import subprocess
import glob
import logging
import logging.handlers
from collections import defaultdict
import re
import textwrap
import datetime
import multiprocessing
from multiprocessing import Manager
from concurrent.futures import wait, as_completed, TimeoutError, CancelledError
import pkg_resources
import time
import signal


# f-strings will produce a 'SyntaxError: invalid syntax' error if not supported by Python version:
f'HybPiper requires Python 3.6 or higher.'

# Import non-standard-library modules:
unsuccessful_imports = []
try:
    from Bio import SeqIO, SeqRecord
except ImportError:
    unsuccessful_imports.append('Bio')

if unsuccessful_imports:
    package_list = '\n'.join(unsuccessful_imports)
    sys.exit(f'The required Python packages are not found:\n\n{package_list}\n\nAre they installed for the Python '
             f'installation used to run HybPiper?')

# Import program modules:
from yang_and_smith import paralogy_subparsers
from yang_and_smith import check_and_batch
from yang_and_smith import align_and_clean


########################################################################################################################
# Define functions
########################################################################################################################

def setup_logger(name, log_file, console_level=logging.INFO, file_level=logging.DEBUG,
                 logger_object_level=logging.DEBUG):
    """
    Function to create a logger instance.

    By default, logs level DEBUG and above to file.
    By default, logs level INFO and above to stderr and file.

    :param str name: name for the logger instance
    :param str log_file: filename for log file
    :param str console_level: logger level for logging to console
    :param str file_level: logger level for logging to file
    :param str logger_object_level: logger level for logger object
    :return: logging.Logger: logger object
    """

    # Get date and time string for log filename:
    date_and_time = datetime.datetime.now().strftime("%Y-%m-%d-%H_%M_%S")

    # Log to file:
    file_handler = logging.FileHandler(f'{log_file}_{date_and_time}.log', mode='w')
    file_handler.setLevel(file_level)
    file_format = logging.Formatter('%(asctime)s - %(filename)s - %(name)s - %(funcName)s - %(levelname)s - %('
                                    'message)s')
    file_handler.setFormatter(file_format)

    # Log to Terminal (stderr):
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


def check_and_batch_main(args):
    """
    Calls the function main() from module check_and_batch

    :param args: argparse namespace with subparser options for function check_and_batch.main()
    :return: None: no return value specified; default is None
    """

    check_and_batch.main(args)


def align_and_clean_main(args):
    """
    Calls the function main() from module align_and_clean

    :param args: argparse namespace with subparser options for function align_and_clean.main()
    :return: None: no return value specified; default is None
    """

    align_and_clean.main(args)


def parse_arguments():
    """
    Creates main parser and add subparsers. Parses command line arguments

    :return argparse.Namespace arguments: arguments for the given command/subcommand
    """

    parser = argparse.ArgumentParser(prog='resolve_paralogs', description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='To view parameters and help for a subcommand, use e.g. "check_and_batch'
                                            '--help"')
    group_1 = parser.add_mutually_exclusive_group(required=False)
    group_1.add_argument('--version', '-v',
                         dest='version',
                         action='version',
                         version='%(prog)s 1.0.0rc build 1',
                         help='Print the resolve_paralogs version number.')

    # Add subparsers:
    subparsers = parser.add_subparsers(title='Subcommands for resolve_paralogs', description='Valid subcommands:')
    parser_check_and_batch = paralogy_subparsers.add_check_and_batch_parser(subparsers)
    parser_align_and_clean = paralogy_subparsers.add_align_and_clean_parser(subparsers)

    # Set functions for subparsers:
    parser_check_and_batch.set_defaults(func=check_and_batch_main)
    parser_align_and_clean.set_defaults(func=align_and_clean_main)

    # Parse and return all arguments:
    arguments = parser.parse_args()

    return arguments


def main():

    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(1)

    # check for external dependencies:


    # Parse arguments for the command/subcommand used:
    args = parse_arguments()

    # Run the function associated with the subcommand:
    args.func(args)


########################################################################################################################
# Run the script
#######################################################################################################################
if __name__ == '__main__':
    main()
    # import pstats
    # import cProfile
    # profiler = cProfile.Profile()
    # profiler.enable()
    # main()
    # profiler.disable()
    # stats = pstats.Stats(profiler).sort_stats('cumtime')
    # stats.print_stats()

################################################## END OF SCRIPT #######################################################