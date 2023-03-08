#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
ParaGone: paralogy resolution pipeline version 0.0.7rc (March 2023)

Adapted from Yang and Smith, Mol Biol Evol. 2014 Nov; 31(11): 3081–3092.

Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

For a list of available subcommands, run:

    paragone --help

"""

import argparse
import sys
import cProfile
import textwrap
import platform

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
from paragone import trim_trees_treeshrink
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

    # Log system details for debugging:
    utils.get_platform_info(logger=logger)

    if platform.system() == 'Darwin':
        utils.check_macos_version(logger=logger)

    # Log ulimit details for debugging:
    utils.get_ulimit_info(logger=logger)

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

    # Log system details for debugging:
    utils.get_platform_info(logger=logger)

    if platform.system() == 'Darwin':
        utils.check_macos_version(logger=logger)

    # Log ulimit details for debugging:
    utils.get_ulimit_info(logger=logger)

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

    # Log system details for debugging:
    utils.get_platform_info(logger=logger)

    if platform.system() == 'Darwin':
        utils.check_macos_version(logger=logger)

    # Log ulimit details for debugging:
    utils.get_ulimit_info(logger=logger)

    # Trim tips from the input tree using TreeShrink:
    trim_trees_treeshrink.main(
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
    logger = utils.setup_logger(__name__, f'{log_directory}/align_selected_and_tree')

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

    # Log system details for debugging:
    utils.get_platform_info(logger=logger)

    if platform.system() == 'Darwin':
        utils.check_macos_version(logger=logger)

    # Log ulimit details for debugging:
    utils.get_ulimit_info(logger=logger)

    align_selected_and_tree.main(
        args,
        report_directory,
        logger=logger)


def prune_paralogs_main(args,
                        log_directory=None,
                        report_directory=None):
    """
    Calls the function main() from module prune_paralogs_rt, prune_paralogs_mo, prune_paralogs_mi

    :param args: argparse namespace with subparser options for function prune_paralogs.main()
    :param str log_directory: path to directory for log files
    :param str report_directory: path to directory for report files
    :return: None: no return value specified; default is None
    """

    # Create a dictionary from the argparse Namespace:
    parameters = vars(args)

    # Create a logger for alignment_to_tree_main:
    logger = utils.setup_logger(__name__, f'{log_directory}/prune_paralogs')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand prune_paralogs was called with these arguments:\n')

    for parameter, value in parameters.items():
        if not parameter == 'func':
            logger.info(f'{" " * 10} {parameter}: {value}')
    logger.info('')

    # Log system details for debugging:
    utils.get_platform_info(logger=logger)

    if platform.system() == 'Darwin':
        utils.check_macos_version(logger=logger)

    # Log ulimit details for debugging:
    utils.get_ulimit_info(logger=logger)

    # Check if at least one pruning algorithm was provided:
    if not args.mo and not args.mi and not args.rt:
        logger.error(f'{"[ERROR]:":10} Please provide at least one of the following parameters: --rt, --mo, --mi\n')
        sys.exit()

    # Run the Monophyletic Outgroups (MO) algorithm:
    if args.mo:
        prune_paralogs_mo.main(
            args,
            report_directory,
            logger=logger)

    # Run the Maximum Inclusion (MI) algorithm:
    if args.mi:
        prune_paralogs_mi.main(
            args,
            report_directory,
            logger=logger)

    # Run the RooTed outgroups (RT) algorithm:
    if args.rt:
        prune_paralogs_rt.main(
            args,
            report_directory,
            logger=logger)


def final_alignments_main(args,
                          log_directory=None,
                          report_directory=None):
    """
    Calls the function main() from module strip_names_and_align

    :param args: argparse namespace with subparser options for function strip_names_and_align.main()
    :param str log_directory: path to directory for log files
    :param str report_directory: path to directory for report files
    :return: None: no return value specified; default is None
    """

    # Create a dictionary from the argparse Namespace:
    parameters = vars(args)

    # Create a logger for alignment_to_tree_main:
    logger = utils.setup_logger(__name__, f'{log_directory}/strip_names_and_align')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand strip_names_and_align was called with these arguments:\n')

    for parameter, value in parameters.items():
        if not parameter == 'func':
            logger.info(f'{" " * 10} {parameter}: {value}')
    logger.info('')

    # Log system details for debugging:
    utils.get_platform_info(logger=logger)

    if platform.system() == 'Darwin':
        utils.check_macos_version(logger=logger)

    # Log ulimit details for debugging:
    utils.get_ulimit_info(logger=logger)

    # Extract fasta for the Monophyletic Outgroups (MO) algorithm:
    if args.mo:
        fasta_from_tree.main(
            args,
            report_directory,
            algorithm_suffix='mo',
            logger=logger)

        strip_names_and_align.main(
            args,
            report_directory,
            selected_alignment_directory='17_selected_sequences_MO',
            logger=logger)

    # Extract fasta for the Maximum Inclusion (MI) algorithm:
    if args.mi:
        fasta_from_tree.main(
            args,
            report_directory,
            algorithm_suffix='mi',
            logger=logger)

        strip_names_and_align.main(
            args,
            report_directory,
            selected_alignment_directory='18_selected_sequences_MI',
            logger=logger)

    # Extract fasta for the RooTed outgroups (RT) algorithm:
    if args.rt:
        fasta_from_tree.main(
            args,
            report_directory,
            algorithm_suffix='rt',
            logger=logger)

        strip_names_and_align.main(
            args,
            report_directory,
            selected_alignment_directory='19_selected_sequences_RT',
            logger=logger)

    # Remove all intermediate files and folders:
    if not args.keep_intermediate_files:
        utils.delete_intermediate_data(logger=logger)


def full_pipeline_main(args,
                       log_directory=None,
                       report_directory=None):
    """
    Runs all steps of the ParaGone pipeline

    :param args: argparse namespace with subparser options for function full_pipeline_main()
    :param str log_directory: path to directory for log files
    :param str report_directory: path to directory for report files
    :return: None: no return value specified; default is None
    """

    # Create a dictionary from the argparse Namespace:
    parameters = vars(args)

    # Create a logger for alignment_to_tree_main:
    logger = utils.setup_logger(__name__, f'{log_directory}/full_pipeline')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand paragone_full_pipeline was called with these arguments:\n')

    for parameter, value in parameters.items():
        if not parameter == 'func':
            logger.info(f'{" " * 10} {parameter}: {value}')
    logger.info('')

    # Log system details for debugging:
    utils.get_platform_info(logger=logger)

    if platform.system() == 'Darwin':
        utils.check_macos_version(logger=logger)

    # Log ulimit details for debugging:
    utils.get_ulimit_info(logger=logger)

    # Check input files:
    check_inputs.main(
        args,
        logger=logger)

    # Align paralogs, trim and clean:
    align_and_clean.main(
        args,
        logger=logger)

    if args.no_trimming and args.no_cleaning:
        args.alignment_directory = '02_alignments'
    elif args.no_cleaning and not args.no_trimming:
        args.alignment_directory = '03_alignments_trimmed'
    elif args.no_trimming and not args.no_cleaning:
        args.alignment_directory = '04_alignments_hmmcleaned'
    else:
        args.alignment_directory = '04_alignments_trimmed_hmmcleaned'

    # Produce trees from alignments:
    alignment_to_tree.main(
        args,
        logger=logger)

    # Trim tips from the input tree using TreeShrink:
    trim_trees_treeshrink.main(
        args,
        report_directory,
        logger=logger)

    # Mask tree tips in trimmed trees:
    args.mask_tips_alignment_directory = args.alignment_directory
    args.mask_tips_alignment_file_suffix = '.fasta'
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
    args.from_cut_deep_paralogs = True
    fasta_from_tree.main(
        args,
        report_directory,
        logger=logger)

    # Align the selected sequences and generate trees:
    args.qc_alignment_directory = args.alignment_directory
    align_selected_and_tree.main(
        args,
        report_directory,
        logger=logger)

    # Run the Monophyletic Outgroups (MO) algorithm:
    if args.mo:
        prune_paralogs_mo.main(
            args,
            report_directory,
            logger=logger)

    # Run the Maximum Inclusion (MI) algorithm:
    if args.mi:
        args.relative_tip_cutoff = args.trim_tips_relative_cutoff
        args.absolute_tip_cutoff = args.trim_tips_absolute_cutoff
        prune_paralogs_mi.main(
            args,
            report_directory,
            logger=logger)

    # Run the RooTed outgroups (RT) algorithm:
    if args.rt:
        prune_paralogs_rt.main(
            args,
            report_directory,
            logger=logger)

    # Extract fasta for the Monophyletic Outgroups (MO) algorithm:
    if args.mo:
        args.from_cut_deep_paralogs = False
        fasta_from_tree.main(
            args,
            report_directory,
            algorithm_suffix='mo',
            logger=logger)

        strip_names_and_align.main(
            args,
            report_directory,
            selected_alignment_directory='17_selected_sequences_MO',
            logger=logger)

    # Extract fasta for the Maximum Inclusion (MI) algorithm:
    if args.mi:
        args.from_cut_deep_paralogs = False
        fasta_from_tree.main(
            args,
            report_directory,
            algorithm_suffix='mi',
            logger=logger)

        strip_names_and_align.main(
            args,
            report_directory,
            selected_alignment_directory='18_selected_sequences_MI',
            logger=logger)

    # Extract fasta for the RooTed outgroups (RT) algorithm:
    if args.rt:
        args.from_cut_deep_paralogs = False
        fasta_from_tree.main(
            args,
            report_directory,
            algorithm_suffix='rt',
            logger=logger)

        strip_names_and_align.main(
            args,
            report_directory,
            selected_alignment_directory='19_selected_sequences_RT',
            logger=logger)

    # Remove all intermediate files and folders:
    if not args.keep_intermediate_files:
        utils.delete_intermediate_data(logger=logger)


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
                         version='%(prog)s 0.0.7rc',
                         help='Print the ParaGone version number.')

    # Add subparsers:
    subparsers = parser.add_subparsers(title='Subcommands for paragone', metavar='', )
    parser_check_and_align = paragone_subparsers.add_check_and_align_parser(subparsers)
    parser_alignment_to_tree = paragone_subparsers.add_alignment_to_tree_parser(subparsers)
    parser_qc_trees_and_extract_fasta = paragone_subparsers.add_qc_trees_and_extract_fasta(subparsers)
    parser_align_selected_and_tree = paragone_subparsers.add_align_selected_and_tree_parser(subparsers)
    parser_prune_paralogs = paragone_subparsers.add_prune_paralogs_parser(subparsers)
    parser_final_alignments = paragone_subparsers.add_final_alignments_parser(subparsers)
    parser_full_pipeline = paragone_subparsers.add_full_pipeline_parser(subparsers)

    # Set functions for subparsers:
    parser_check_and_align.set_defaults(func=check_and_align_main)
    parser_alignment_to_tree.set_defaults(func=alignment_to_tree_main)
    parser_qc_trees_and_extract_fasta.set_defaults(func=qc_trees_and_extract_fasta_main)
    parser_align_selected_and_tree.set_defaults(func=align_selected_and_tree_main)
    parser_prune_paralogs.set_defaults(func=prune_paralogs_main)
    parser_final_alignments.set_defaults(func=final_alignments_main)
    parser_full_pipeline.set_defaults(func=full_pipeline_main)

    # Parse and return all arguments:
    arguments = parser.parse_args()

    return arguments


def main():

    title = textwrap.dedent(

        fr"""

     _____      ATCTATCTATAC.......       ___        
    |  _  \                             /  ___\            
    | |_| |   _____   _____     _____  / /   __    _____    _   __     _____              
    |  ___/ /  _   | |  _  \  /  _   | | | |_  \  /  _  \  | |/   \   |  _  |        
    | |    |  |_|  | | |  -- |  |_|  | | |___| | |  |_|  | |    \  \  |  __/              
    |_|     \ __/ _| |_|      \ __/ _|  \ ____ /  \_____/  |_|   |__| |_____|      

                                             .......ATCGACTGCACGTGACTCG        

        """)

    sys.stderr.write(title)
    sys.stderr.flush()

    if len(sys.argv) == 1:
        sys.stderr.write(__doc__)
        sys.exit(1)

    # Parse arguments for the command/subcommand used:
    args = parse_arguments()

    # Create a directory for logs for each step of the pipeline:
    utils.createfolder('00_logs_and_reports')
    log_directory = utils.createfolder('00_logs_and_reports/logs')
    report_directory = utils.createfolder('00_logs_and_reports/reports')

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
