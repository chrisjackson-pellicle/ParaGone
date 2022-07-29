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
try:
    from ete3 import Tree
except ImportError:
    unsuccessful_imports.append('ete3')


if unsuccessful_imports:
    package_list = '\n'.join(unsuccessful_imports)
    sys.exit(f'The required Python packages are not found:\n\n{package_list}\n\nAre they installed for the Python '
             f'installation used to run resolve_paralogs?')

# Import program modules:
from yang_and_smith import paralogy_subparsers
from yang_and_smith import check_and_batch
from yang_and_smith import align_and_clean
from yang_and_smith import alignment_to_tree
from yang_and_smith import trim_tree_tips
from yang_and_smith import mask_tree_tips
from yang_and_smith import cut_deep_paralogs
from yang_and_smith import fasta_from_tree
from yang_and_smith import newick3
from yang_and_smith import phylo3
from yang_and_smith import tree_utils
from yang_and_smith import seq
from yang_and_smith import utils


########################################################################################################################
# Define functions
########################################################################################################################


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


def alignment_to_tree_main(args):
    """
    Calls the function main() from module alignment_to_tree

    :param args: argparse namespace with subparser options for function alignment_to_tree.main()
    :return: None: no return value specified; default is None
    """

    alignment_to_tree.main(args)


def trim_tree_tips_main(args):
    """
    Calls the function main() from module trim_tree_tips

    :param args: argparse namespace with subparser options for function trim_tree_tips.main()
    :return: None: no return value specified; default is None
    """

    trim_tree_tips.main(args)


def mask_tree_tips_main(args):
    """
    Calls the function main() from module mask_tree_tips

    :param args: argparse namespace with subparser options for function mask_tree_tips.main()
    :return: None: no return value specified; default is None
    """

    mask_tree_tips.main(args)


def cut_deep_paralogs_main(args):
    """
    Calls the function main() from module cut_deep_paralogs

    :param args: argparse namespace with subparser options for function cut_deep_paralogs.main()
    :return: None: no return value specified; default is None
    """

    cut_deep_paralogs.main(args)


def fasta_from_tree_main(args):
    """
    Calls the function main() from module fasta_from_tree

    :param args: argparse namespace with subparser options for function fasta_from_tree.main()
    :return: None: no return value specified; default is None
    """

    fasta_from_tree.main(args)


def parse_arguments():
    """
    Creates main parser and add subparsers. Parses command line arguments

    :return argparse.Namespace arguments: arguments for the given command/subcommand
    """

    parser = argparse.ArgumentParser(prog='resolve_paralogs', description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='To view parameters and help for a subcommand, use e.g. "check_and_batch '
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
    parser_alignment_to_tree = paralogy_subparsers.add_alignment_to_tree_parser(subparsers)
    parser_trim_tree_tips = paralogy_subparsers.add_trim_tree_tips_parser(subparsers)
    parser_mask_tree_tips = paralogy_subparsers.add_mask_tree_tips_parser(subparsers)
    parser_cut_deep_paralogs = paralogy_subparsers.add_cut_deep_paralogs_parser(subparsers)
    parser_fasta_from_tree = paralogy_subparsers.add_fasta_from_tree_parser(subparsers)

    # Set functions for subparsers:
    parser_check_and_batch.set_defaults(func=check_and_batch_main)
    parser_align_and_clean.set_defaults(func=align_and_clean_main)
    parser_alignment_to_tree.set_defaults(func=alignment_to_tree_main)
    parser_trim_tree_tips.set_defaults(func=trim_tree_tips_main)
    parser_mask_tree_tips.set_defaults(func=mask_tree_tips_main)
    parser_cut_deep_paralogs.set_defaults(func=cut_deep_paralogs_main)
    parser_fasta_from_tree.set_defaults(func=fasta_from_tree_main)

    # Parse and return all arguments:
    arguments = parser.parse_args()

    return arguments


def main():

    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(1)

    # Create a directory for logs for each step of the pipeline:
    utils.createfolder('logs_resolve_paralogs')

    # Initialise logger:
    # logger = utils.setup_logger(__name__, 'logs_resolve_paralogs/resolve_paralogs')

    # check for external dependencies:
    # if utils.check_dependencies(logger=logger):
    #     logger.info(f'{"[INFO]:":10} All external dependencies found!')
    # else:
    #     logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
    #     sys.exit(1)

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