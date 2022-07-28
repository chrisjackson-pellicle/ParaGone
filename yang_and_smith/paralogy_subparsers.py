#!/usr/bin/env python

"""
Contains argument subparsers
"""

import textwrap
import logging
import sys


def add_check_and_batch_parser(subparsers):
    """
    Parser for check_and_batch

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_check_and_batch = subparsers.add_parser('check_and_batch',
                                                   help='Check input file, outgroup coverage, create batches')
    parser_check_and_batch.add_argument('gene_fasta_directory',
                                        type=str,
                                        help="Directory contains fasta files including paralogs")
    parser_check_and_batch.add_argument('--external_outgroups_file',
                                        type=str,
                                        default=None,
                                        help='File in fasta format with additional outgroup sequences to add to each '
                                             'gene')
    parser_check_and_batch.add_argument('--external_outgroup',
                                        action='append',
                                        type=str,
                                        dest='external_outgroups',
                                        default=None,
                                        help='If a taxon name is provided, only use these sequences from '
                                             'the user-provided external_outgroups_file. Note that this parameter can '
                                             'be specified one ore more times.')
    parser_check_and_batch.add_argument('--internal_outgroup',
                                        action='append',
                                        type=str,
                                        dest='internal_outgroups',
                                        default=None,
                                        help='Taxon name to use as an internal outgroup (i.e. present in input '
                                             'paralog files). Note that this parameter can be specified one ore more '
                                             'times.')
    parser_check_and_batch.add_argument('--batch_size',
                                        type=int,
                                        default=20,
                                        help='Number of fasta files for each batch of input paralog fasta files. '
                                             'Default is: %(default)s')

    return parser_check_and_batch


def add_align_and_clean_parser(subparsers):
    """
    Parser for align_and_clean

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_align_and_clean = subparsers.add_parser('align_and_clean',
                                                   help='Align and QC input files')
    parser_align_and_clean.add_argument('gene_fasta_directory',
                                        type=str,
                                        help='directory containing fasta files (with sanitised gene names) including '
                                             'paralogs')
    parser_align_and_clean.add_argument('--pool',
                                        type=int,
                                        default=1,
                                        help='Number of alignments to run concurrently. Default is: %(default)s')
    parser_align_and_clean.add_argument('--threads',
                                        type=int,
                                        default=1,
                                        help='Number of threads to use for each concurrent alignment. Default '
                                             'is: %(default)s')
    parser_align_and_clean.add_argument('--no_stitched_contigs',
                                        action='store_true',
                                        default=False,
                                        help='If specified, realign mafft alignments with clustal omega. Default is: '
                                             '%(default)s')
    parser_align_and_clean.add_argument('--use_muscle',
                                        action='store_true',
                                        default=False,
                                        help='If specified, use muscle rather than mafft for initial alignments. '
                                             'Default is: %(default)s')
    parser_align_and_clean.add_argument('--mafft_algorithm',
                                        default='auto',
                                        help='Algorithm to use for mafft alignments. Default is: %(default)s')

    return parser_align_and_clean


def add_alignment_to_tree_parser(subparsers):
    """
    Parser for alignment_to_tree

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_alignment_to_tree = subparsers.add_parser('alignment_to_tree',
                                                     help='Takes folder of alignments and generates phylogenetic trees')
    parser_alignment_to_tree.add_argument('alignment_directory',
                                          type=str,
                                          help='directory containing fasta alignment files')
    parser_alignment_to_tree.add_argument('--pool',
                                          type=int,
                                          default=1,
                                          help='Number of trees to run concurrently. Default is: %(default)s')
    parser_alignment_to_tree.add_argument('--threads',
                                          type=int,
                                          default=1,
                                          help='Number of threads to use for each concurrent tree. Default '
                                               'is: %(default)s')
    parser_alignment_to_tree.add_argument('--generate_bootstraps',
                                          action='store_true',
                                          default=False,
                                          help='Create bootstraps for trees using UFBoot. Default is: '
                                               '%(default)s')
    parser_alignment_to_tree.add_argument('--use_fasttree',
                                          action='store_true',
                                          default=False,
                                          help='Use FastTree instead of IQTREE. Default is: %(default)s')

    return parser_alignment_to_tree


def add_trim_tree_tips_parser(subparsers):
    """
    Parser for trim_tree_tips

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_trim_tree_tips = subparsers.add_parser('trim_tree_tips',
                                                  help='Remove tree tips on long branches')
    parser_trim_tree_tips.add_argument('treefile_directory',
                                       type=str,
                                       help='directory containing tree newick files')
    parser_trim_tree_tips.add_argument('--tree_file_suffix',
                                       type=str,
                                       default='.treefile',
                                       help='Suffix for newick tree files. Default is: %(default)s')
    parser_trim_tree_tips.add_argument('--relative_cutoff',
                                       type=float,
                                       default=0.2,
                                       help='Relative cutoff for removing tree tips. Default is: %(default)s')
    parser_trim_tree_tips.add_argument('--absolute_cutoff',
                                       type=float,
                                       default=0.4,
                                       help='Absolute cutoff for removing tree tips. Default is: %(default)s')
    parser_trim_tree_tips.add_argument('--output_folder',
                                       default='tree_files_trimmed',
                                       help='Directory name for output trees. Default is: %(default)s')

    return parser_trim_tree_tips



# def add_mask_tree_tips_parser(subparsers):
#     """
#     Parser for mask_tree_tips
#
#     :param argparse._SubParsersAction subparsers:
#     :return None: no return value specified; default is None
#     """
#
#     parser_mask_tree_tips = subparsers.add_parser('mask_tree_tips',
#                                                    help='Takes folder of alignments and generates phylogenetic trees')
#     parser_mask_tree_tips.add_argument('alignment_directory',
#                                           type=str,
#                                           help='directory containing fasta alignment files')
#     parser_mask_tree_tips.add_argument('--pool',
#                                           type=int,
#                                           default=1,
#                                           help='Number of trees to run concurrently. Default is: %(default)s')
#     parser_mask_tree_tips.add_argument('--threads',
#                                           type=int,
#                                           default=1,
#                                           help='Number of threads to use for each concurrent tree. Default '
#                                                'is: %(default)s')
#     parser_mask_tree_tips.add_argument('--generate_bootstraps',
#                                           action='store_true',
#                                           default=False,
#                                           help='Create bootstraps for trees using UFBoot. Default is: '
#                                                '%(default)s')
#     parser_mask_tree_tips.add_argument('--use_fasttree',
#                                           action='store_true',
#                                           default=False,
#                                           help='Use FastTree instead of IQTREE. Default is: %(default)s')
#
#     return parser_mask_tree_tips
