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


def check_and_align_main(args):
    """
    Calls the function main() from modules check_inputs followed by align_and_clean

    :param args: argparse namespace with subparser options for function check_and_align.main()
    :return: None: no return value specified; default is None
    """

    check_inputs.main(args)
    align_and_clean.main(args)


def alignment_to_tree_main(args):
    """
    Calls the function main() from module alignment_to_tree

    :param args: argparse namespace with subparser options for function alignment_to_tree.main()
    :return: None: no return value specified; default is None
    """

    alignment_to_tree.main(args)


def collate_alignments_and_trees_main(args):
    """
    Calls the function main() from module collate_alignments_and_trees

    :param args: argparse namespace with subparser options for function collate_alignments_and_trees.main()
    :return: None: no return value specified; default is None
    """

    collate_alignments_and_trees.main(args)


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


def align_selected_and_tree_main(args):
    """
    Calls the function main() from module align_selected_and_tree

    :param args: argparse namespace with subparser options for function align_selected_and_tree.main()
    :return: None: no return value specified; default is None
    """

    align_selected_and_tree.main(args)


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
    # parser_alignment_to_tree = paragone_subparsers.add_alignment_to_tree_parser(subparsers)
    # parser_collate_alignments_and_trees = paragone_subparsers.add_collate_alignments_and_trees_parser(subparsers)
    # parser_trim_tree_tips = paragone_subparsers.add_trim_tree_tips_parser(subparsers)
    # parser_mask_tree_tips = paragone_subparsers.add_mask_tree_tips_parser(subparsers)
    # parser_cut_deep_paralogs = paragone_subparsers.add_cut_deep_paralogs_parser(subparsers)
    # parser_fasta_from_tree = paragone_subparsers.add_fasta_from_tree_parser(subparsers)
    # parser_align_selected_and_tree = paragone_subparsers.add_align_selected_and_tree_parser(subparsers)
    # parser_prune_paralogs_mo = paragone_subparsers.add_prune_paralogs_mo_parser(subparsers)
    # parser_prune_paralogs_rt = paragone_subparsers.add_prune_paralogs_rt_parser(subparsers)
    # parser_prune_paralogs_mi = paragone_subparsers.add_prune_paralogs_mi_parser(subparsers)
    # parser_strip_names_and_align = paragone_subparsers.add_strip_names_and_align_parser(subparsers)

    # Set functions for subparsers:
    # parser_check_and_batch.set_defaults(func=check_and_batch_main)
    # parser_align_and_clean.set_defaults(func=align_and_clean_main)
    parser_check_and_align.set_defaults(func=check_and_align_main)
    # parser_alignment_to_tree.set_defaults(func=alignment_to_tree_main)
    # parser_collate_alignments_and_trees.set_defaults(func=collate_alignments_and_trees_main)
    # parser_trim_tree_tips.set_defaults(func=trim_tree_tips_main)
    # parser_mask_tree_tips.set_defaults(func=mask_tree_tips_main)
    # parser_cut_deep_paralogs.set_defaults(func=cut_deep_paralogs_main)
    # parser_fasta_from_tree.set_defaults(func=fasta_from_tree_main)
    # parser_align_selected_and_tree.set_defaults(func=align_selected_and_tree_main)
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
    utils.createfolder('00_logs_and_reports_resolve_paralogs/logs')
    utils.createfolder('00_logs_and_reports_resolve_paralogs/reports')

    # Parse arguments for the command/subcommand used:
    args = parse_arguments()

    # Run the function associated with the subcommand, with or without cProfile:
    if args.run_profiler:
        profiler = cProfile.Profile()
        profiler.enable()
        args.func(args)
        profiler.disable()
        csv = utils.cprofile_to_csv(profiler)

        with open(f'{sys.argv[1]}_cprofile.csv', 'w+') as cprofile_handle:
            cprofile_handle.write(csv)
    else:
        args.func(args)


########################################################################################################################
# Run the script
#######################################################################################################################
if __name__ == '__main__':
    main()

################################################## END OF SCRIPT #######################################################