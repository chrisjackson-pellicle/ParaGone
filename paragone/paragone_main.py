#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
ParaGone: paralogy resolution pipeline version 1.0.0 release candidate (November 2022)

Adapted from Yang and Smith, Mol Biol Evol. 2014 Nov; 31(11): 3081–3092.

Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

For a list of available subcommands, run:

    paragone --help

"""

import argparse
import sys
import cProfile
import logging
import textwrap

# f-strings will produce a 'SyntaxError: invalid syntax' error if not supported by Python version:
f'ParaGone requires Python 3.6 or higher.'

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
             f'installation used to run paragone?')

# Import program modules:
from paragone import paragone_subparsers
from paragone import check_inputs
from paragone import align_and_clean
from paragone import alignment_to_tree
from paragone import collate_alignments_and_trees
from paragone import trim_tree_tips
from paragone import mask_tree_tips
from paragone import cut_deep_paralogs
from paragone import fasta_from_tree
from paragone import align_selected_and_tree
from paragone import prune_paralogs_mo
from paragone import prune_paralogs_rt
from paragone import prune_paralogs_mi
from paragone import strip_names_and_align
from paragone import utils

########################################################################################################################
# Define functions
########################################################################################################################


def check_and_align_main(args,
                         log_directory=None,
                         report_directory=None):
    """
    Calls the function main() from modules: check_inputs -> align_and_clean

    :param args: argparse namespace with subparser options for function check_and_align.main()
    :param str log_directory: path to directory for log files
    :param str report_directory: path to directory for report files
    :return: None: no return value specified; default is None
    """

    # Create a dictionary from the argparse Namespace:
    parameters = vars(args)

    # Create a logger for check_and_align_main:
    logger = utils.setup_logger(__name__, f'{log_directory}/check_and_align')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand check_and_align was called with these arguments:\n')

    for parameter, value in parameters.items():
        if not parameter == 'func':
            logger.info(f'{" "* 10} {parameter}: {value}')
    logger.info('')

    # Check input files:
    check_inputs.main(
        args,
        logger=logger)

    # Align paralogs, trim and clean:
    align_and_clean.main(
        args,
        logger=logger)


def alignment_to_tree_main(args,
                           log_directory=None,
                           report_directory=None):
    """
    Calls the function main() from module alignment_to_tree

    :param args: argparse namespace with subparser options for function alignment_to_tree.main()
    :param str log_directory: path to directory for log files
    :param str report_directory: path to directory for report files
    :return: None: no return value specified; default is None
    """

    # Create a dictionary from the argparse Namespace:
    parameters = vars(args)

    # Create a logger for alignment_to_tree_main:
    logger = utils.setup_logger(__name__, f'{log_directory}/alignment_to_tree')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand alignment_to_tree was called with these arguments:\n')

    for parameter, value in parameters.items():
        if not parameter == 'func':
            logger.info(f'{" " * 10} {parameter}: {value}')
    logger.info('')

    # Produce trees from alignments:
    alignment_to_tree.main(
        args,
        logger=logger)


def qc_trees_and_extract_fasta_main(args,
                                    log_directory=None,
                                    report_directory=None):
    """
    Calls the function main() from modules: trim_tree_tips -> mask_tree_tips -> cut_deep_paralogs -> fasta_from_tree

    :param args: argparse namespace with subparser options for function qc_trees_and_fasta_main.main()
    :param str log_directory: path to directory for log files
    :param str report_directory: path to directory for report files
    :return: None: no return value specified; default is None
    """

    # Create a dictionary from the argparse Namespace:
    parameters = vars(args)

    # Create a logger for qc_trees_and_extract_fasta_main:
    logger = utils.setup_logger(__name__, f'{log_directory}/qc_trees_and_extract_fasta')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand qc_trees_and_extract_fasta was called with these arguments:\n')

    for parameter, value in parameters.items():
        if parameter not in ['func', 'from_cut_deep_paralogs']:
            logger.info(f'{" " * 10} {parameter}: {value}')
    logger.info('')

    # Trim tips from the input tree if above threshold lengths:
    trim_tree_tips.main(
        args,
        report_directory,
        logger=logger)

    # Mask tree tips in trimmed trees:
    mask_tree_tips.main(
        args,
        report_directory,
        logger=logger)

    # Cut tree at points of putative deep paralogs:
    cut_deep_paralogs.main(
        args,
        report_directory,
        logger=logger)

    # Extract fasta sequences corresponding to the QC'd trees:
    fasta_from_tree.main(
        args,
        report_directory,
        logger=logger)


def align_selected_and_tree_main(args,
                                 log_directory=None,
                                 report_directory=None):
    """
    Calls the function main() from module align_selected_and_tree

    :param args: argparse namespace with subparser options for function align_selected_and_tree.main()
    :param str log_directory: path to directory for log files
    :param str report_directory: path to directory for report files
    :return: None: no return value specified; default is None
    """

    # Create a dictionary from the argparse Namespace:
    parameters = vars(args)

    # Create a logger for qc_trees_and_extract_fasta_main:
    logger = utils.setup_logger(__name__, f'{log_directory}/qc_trees_and_extract_fasta')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand qc_trees_and_extract_fasta was called with these arguments:\n')

    for parameter, value in parameters.items():
        if parameter not in ['func', 'from_cut_deep_paralogs']:
            logger.info(f'{" " * 10} {parameter}: {value}')
    logger.info('')

    align_selected_and_tree.main(
        args,
        report_directory,
        logger=logger)


def prune_paralogs_mo_main(args):
    """
    Calls the function main() from module prune_paralogs_mo

    :param args: argparse namespace with subparser options for function prune_paralogs_mo.main()
    :return: None: no return value specified; default is None
    """

    prune_paralogs_mo.main(args)


def prune_paralogs_rt_main(args):
    """
    Calls the function main() from module prune_paralogs_rt

    :param args: argparse namespace with subparser options for function prune_paralogs_rt.main()
    :return: None: no return value specified; default is None
    """

    prune_paralogs_rt.main(args)


def prune_paralogs_mi_main(args):
    """
    Calls the function main() from module prune_paralogs_mi

    :param args: argparse namespace with subparser options for function prune_paralogs_mi.main()
    :return: None: no return value specified; default is None
    """

    prune_paralogs_mi.main(args)


def strip_names_and_align_main(args):
    """
    Calls the function main() from module strip_names_and_align

    :param args: argparse namespace with subparser options for function strip_names_and_align.main()
    :return: None: no return value specified; default is None
    """

    strip_names_and_align.main(args)


def parse_arguments():
    """
    Creates main parser and add subparsers. Parses command line arguments

    :return argparse.Namespace arguments: arguments for the given command/subcommand
    """

    parser = argparse.ArgumentParser(prog='paragone', description=__doc__,
                                     formatter_class=argparse.RawTextHelpFormatter,
                                     epilog='To view parameters and help for a subcommand, use e.g. "check_and_align '
                                            '--help"')
    group_1 = parser.add_mutually_exclusive_group(required=False)
    group_1.add_argument('--version', '-v',
                         dest='version',
                         action='version',
                         version='%(prog)s 1.0.1rc',
                         help='Print the paragone version number.')

    # Add subparsers:
    subparsers = parser.add_subparsers(title='Subcommands for paragone', description='Valid subcommands:')
    parser_check_and_align = paragone_subparsers.add_check_and_align_parser(subparsers)
    parser_alignment_to_tree = paragone_subparsers.add_alignment_to_tree_parser(subparsers)
    parser_qc_trees_and_extract_fasta = paragone_subparsers.add_qc_trees_and_extract_fasta(subparsers)
    parser_align_selected_and_tree = paragone_subparsers.add_align_selected_and_tree_parser(subparsers)
    # parser_prune_paralogs_mo = paragone_subparsers.add_prune_paralogs_mo_parser(subparsers)
    # parser_prune_paralogs_rt = paragone_subparsers.add_prune_paralogs_rt_parser(subparsers)
    # parser_prune_paralogs_mi = paragone_subparsers.add_prune_paralogs_mi_parser(subparsers)
    # parser_strip_names_and_align = paragone_subparsers.add_strip_names_and_align_parser(subparsers)

    # Set functions for subparsers:
    parser_check_and_align.set_defaults(func=check_and_align_main)
    parser_alignment_to_tree.set_defaults(func=alignment_to_tree_main)
    parser_qc_trees_and_extract_fasta.set_defaults(func=qc_trees_and_extract_fasta_main)
    parser_align_selected_and_tree.set_defaults(func=align_selected_and_tree_main)
    # parser_prune_paralogs_mo.set_defaults(func=prune_paralogs_mo_main)
    # parser_prune_paralogs_rt.set_defaults(func=prune_paralogs_rt_main)
    # parser_prune_paralogs_mi.set_defaults(func=prune_paralogs_mi_main)
    # parser_strip_names_and_align.set_defaults(func=strip_names_and_align_main)

    # Parse and return all arguments:
    arguments = parser.parse_args()

    return arguments


def main():

    if len(sys.argv) == 1:
        print(__doc__)
        sys.exit(1)

    # Create a directory for logs for each step of the pipeline:
    utils.createfolder('00_logs_and_reports_resolve_paralogs')
    log_directory = utils.createfolder('00_logs_and_reports_resolve_paralogs/logs')
    report_directory = utils.createfolder('00_logs_and_reports_resolve_paralogs/reports')

    # Parse arguments for the command/subcommand used:
    args = parse_arguments()

    # # check for external dependencies:
    # if utils.check_dependencies(logger=logger):
    #     logger.info(f'{"[INFO]:":10} All external dependencies found!')
    # else:
    #     logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
    #     sys.exit(1)

    # Run the function associated with the subcommand, with or without cProfile:
    if args.run_profiler:
        profiler = cProfile.Profile()
        profiler.enable()

        args.func(args,
                  log_directory=log_directory,
                  report_directory=report_directory)

        profiler.disable()
        csv = utils.cprofile_to_csv(profiler)

        with open(f'{sys.argv[1]}_cprofile.csv', 'w+') as cprofile_handle:
            cprofile_handle.write(csv)
    else:
        args.func(args,
                  log_directory=log_directory,
                  report_directory=report_directory)


########################################################################################################################
# Run the script
#######################################################################################################################
if __name__ == '__main__':
    main()

################################################## END OF SCRIPT #######################################################