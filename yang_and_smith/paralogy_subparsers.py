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
                                             'paralog files). Note that this parameter can be specified one or more '
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
    # parser_trim_tree_tips.add_argument('--output_folder',
    #                                    default='tree_files_trimmed',
    #                                    help='Directory name for output trees. Default is: %(default)s')

    return parser_trim_tree_tips


def add_mask_tree_tips_parser(subparsers):
    """
    Parser for mask_tree_tips

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_mask_tree_tips = subparsers.add_parser('mask_tree_tips',
                                                  help='Remove monophyletic tree tips from same taxon')
    parser_mask_tree_tips.add_argument('treefile_directory',
                                       type=str,
                                       help='directory containing tree newick files')
    parser_mask_tree_tips.add_argument('--tree_file_suffix',
                                       type=str,
                                       default='.tt',
                                       help='Suffix for newick tree files. Default is: %(default)s')
    parser_mask_tree_tips.add_argument('alignment_directory',
                                       type=str,
                                       help='directory containing fasta alignment files')
    parser_mask_tree_tips.add_argument('--alignment_file_suffix',
                                       type=str,
                                       default='.fasta',
                                       help='Suffix for alignment files. Default is: %(default)s')
    parser_mask_tree_tips.add_argument('--remove_paraphyletic_tips',
                                       action='store_true',
                                       default=False,
                                       help='Remove paraphyletic tree tips. Default is: %(default)s')
    parser_mask_tree_tips.add_argument('--output_folder',
                                       default='tree_files_trimmed_and_masked',
                                       help='Directory name for output trees. Default is: %(default)s')

    return parser_mask_tree_tips


def add_cut_deep_paralogs_parser(subparsers):
    """
    Parser for cut_deep_paralogs

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_cut_deep_paralogs = subparsers.add_parser('cut_deep_paralogs',
                                                     help='Cut long internal branches to remove deep paralogs')
    parser_cut_deep_paralogs.add_argument('treefile_directory',
                                          type=str,
                                          help='directory containing tree newick files')
    parser_cut_deep_paralogs.add_argument('--tree_file_suffix',
                                          type=str,
                                          default='.mm',
                                          help='Suffix for newick tree files. Default is: %(default)s')
    parser_cut_deep_paralogs.add_argument('--internal_branch_length_cutoff',
                                          type=float,
                                          default=0.3,
                                          help='Internal branch length cutoff cutting tree. Default is: %('
                                               'default)s')
    parser_cut_deep_paralogs.add_argument('--minimum_number_taxa',
                                          type=int,
                                          default=4,
                                          help='Minimum number of taxa in tree for tree to be retained. Default is: %('
                                               'default)s')
    parser_cut_deep_paralogs.add_argument('--output_folder',
                                          default='tree_files_trimmed_and_masked_and_cut',
                                          help='Directory name for output trees. Default is: %(default)s')

    return parser_cut_deep_paralogs


def add_fasta_from_tree_parser(subparsers):
    """
    Parser for fasta_from_tree

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_fasta_from_tree = subparsers.add_parser('fasta_from_tree',
                                                   help='Filters an alignment for names in tree')
    parser_fasta_from_tree.add_argument('treefile_directory',
                                        type=str,
                                        help='directory containing tree newick files')
    parser_fasta_from_tree.add_argument('alignment_directory',
                                        type=str,
                                        help='directory containing fasta alignment files')
    parser_fasta_from_tree.add_argument('--tree_file_suffix',
                                        type=str,
                                        default='.subtree',
                                        help='Suffix for newick tree files. Default is: %(default)s')
    parser_fasta_from_tree.add_argument('--from_cut_deep_paralogs',
                                        action='store_true',
                                        default=False,
                                        help='If set, trees are from QC step "cut_deep_paralogs". Default is: %('
                                             'default)s')
    parser_fasta_from_tree.add_argument('--batch_size',
                                        type=int,
                                        default=20,
                                        help='Number of fasta files in each batch, from input paralog fasta files. '
                                             'Default is: %(default)s')

    return parser_fasta_from_tree


def add_align_selected_and_tree_parser(subparsers):
    """
    Parser for align_selected_and_tree

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_align_selected_and_tree = subparsers.add_parser('align_selected_and_tree',
                                                           help='Aligns selected fasta seqs for each subtree, '
                                                                'and generates a new tree')
    parser_align_selected_and_tree.add_argument('selected_alignment_directory',
                                                type=str,
                                                help='directory containing selected alignment files corresponding to '
                                                     'subtrees')
    parser_align_selected_and_tree.add_argument('hmmcleaned_alignment_directory',
                                                type=str,
                                                help='directory containing hmmcleaned original fasta alignment files')
    parser_align_selected_and_tree.add_argument('--pool',
                                                type=int,
                                                default=1,
                                                help='Number of alignments/trees to run concurrently. Default is: %('
                                                     'default)s')
    parser_align_selected_and_tree.add_argument('--threads',
                                                type=int,
                                                default=1,
                                                help='Number of threads to use for each concurrent alignment/tree. '
                                                     'Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--no_stitched_contigs',
                                                action='store_true',
                                                default=False,
                                                help='If specified, realign mafft alignments with clustal omega. '
                                                     'Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--use_muscle',
                                                action='store_true',
                                                default=False,
                                                help='If specified, use muscle rather than mafft for initial '
                                                     'alignments. Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--mafft_algorithm',
                                                default='auto',
                                                help='Algorithm to use for mafft alignments. Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--external_outgroups_file',
                                                type=str,
                                                default=None,
                                                help='file in fasta format with additional outgroup sequences to add '
                                                     'to each gene')
    parser_align_selected_and_tree.add_argument('--external_outgroup',
                                                action='append',
                                                type=str,
                                                dest='external_outgroups',
                                                default=None,
                                                help='If one or more taxon names are provided, only use these '
                                                     'sequences from the user-provided external_outgroups_file')
    parser_align_selected_and_tree.add_argument('--internal_outgroup',
                                                action='append',
                                                type=str,
                                                dest='internal_outgroups',
                                                default=None,
                                                help='Taxon name to use as an internal outgroup (i.e. present in input '
                                                     'paralog files). Note that this parameter can be specified one '
                                                     'or more times.')
    parser_align_selected_and_tree.add_argument('--generate_bootstraps',
                                                action='store_true',
                                                default=False,
                                                help='Create bootstraps for trees using UFBoot. Default is: %('
                                                     'default)s')
    parser_align_selected_and_tree.add_argument('--use_fasttree',
                                                action='store_true',
                                                default=False,
                                                help='Use FastTree instead of IQTREE. Default is: %(default)s')

    return parser_align_selected_and_tree


def add_prune_paralogs_mo_parser(subparsers):
    """
    Parser for fasta_from_tree

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_prune_paralogs_mo = subparsers.add_parser('prune_paralogs_mo',
                                                     help='Prune paralogs from tree using Monophyletic Outgroups (MO) '
                                                     'algorithm')
    parser_prune_paralogs_mo.add_argument('treefile_directory',
                                          type=str,
                                          help='directory containing tree newick files')
    parser_prune_paralogs_mo.add_argument('in_and_outgroup_list',
                                          type=str,
                                          help='Text file listing in- and out-group taxa')
    parser_prune_paralogs_mo.add_argument('--tree_file_suffix',
                                          type=str,
                                          default='.treefile',
                                          help='Suffix for newick tree files. Default is: %(default)s')
    parser_prune_paralogs_mo.add_argument('--minimum_taxa',
                                          type=int,
                                          default=2,
                                          help='Minimum number of taxa required. Default is: %(default)s')
    parser_prune_paralogs_mo.add_argument('--ignore_1to1_orthologs',
                                          action='store_true',
                                          default=False,
                                          help='Output 1to1 orthologs, i.e. trees with no paralogs. Default is: %('
                                               'default)s')

    return parser_prune_paralogs_mo


def add_prune_paralogs_rt_parser(subparsers):
    """
    Parser for fasta_from_tree

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_prune_paralogs_rt = subparsers.add_parser('prune_paralogs_rt',
                                                     help='Prune paralogs from tree using RooTed ingroups (RT) '
                                                     'algorithm')
    parser_prune_paralogs_rt.add_argument('treefile_directory',
                                          type=str,
                                          help='directory containing tree newick files')
    parser_prune_paralogs_rt.add_argument('in_and_outgroup_list',
                                          type=str,
                                          help='Text file listing in- and out-group taxa')
    parser_prune_paralogs_rt.add_argument('--tree_file_suffix',
                                          type=str,
                                          default='.treefile',
                                          help='Suffix for newick tree files. Default is: %(default)s')
    parser_prune_paralogs_rt.add_argument('--minimum_taxa',
                                          type=int,
                                          default=2,
                                          help='Minimum number of taxa required. Default is: %(default)s')

    return parser_prune_paralogs_rt


def add_prune_paralogs_mi_parser(subparsers):
    """
    Parser for fasta_from_tree

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_prune_paralogs_mi = subparsers.add_parser('prune_paralogs_mi',
                                                     help='Prune paralogs from tree using Maximum Inclusion (MI) '
                                                     'algorithm')
    parser_prune_paralogs_mi.add_argument('treefile_directory',
                                          type=str,
                                          help='directory containing tree newick files')
    parser_prune_paralogs_mi.add_argument('in_and_outgroup_list',
                                          type=str,
                                          help='Text file listing in- and out-group taxa')
    parser_prune_paralogs_mi.add_argument('--tree_file_suffix',
                                          type=str,
                                          default='.treefile',
                                          help='Suffix for newick tree files. Default is: %(default)s')
    parser_prune_paralogs_mi.add_argument('--relative_tip_cutoff',
                                          type=float,
                                          default=0.2,
                                          help='Relative tip cut-off threshold. Default is: %(default)s')
    parser_prune_paralogs_mi.add_argument('--absolute_tip_cutoff',
                                          type=float,
                                          default=0.4,
                                          help='Absolute tip cut-off threshold. Default is: %(default)s')
    parser_prune_paralogs_mi.add_argument('--minimum_taxa',
                                          type=int,
                                          default=2,
                                          help='Minimum number of taxa required. Default is: %(default)s')
    parser_prune_paralogs_mi.add_argument('--ignore_1to1_orthologs',
                                          action='store_true',
                                          default=False,
                                          help='Output 1to1 orthologs, i.e. trees with no paralogs. Default is: %('
                                               'default)s')

    return parser_prune_paralogs_mi


def add_strip_names_and_align_parser(subparsers):
    """
    Parser for strip_names_and_align

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_align_selected_and_tree = subparsers.add_parser('strip_names_and_align',
                                                           help='XXX')
    parser_align_selected_and_tree.add_argument('selected_alignment_directory',
                                                type=str,
                                                help='directory containing selected alignment files corresponding to '
                                                     'pruned trees from one of MO, RT, MI algorithms.')
    parser_align_selected_and_tree.add_argument('--pool',
                                                type=int,
                                                default=1,
                                                help='Number of alignments to run concurrently. Default is: %('
                                                     'default)s')
    parser_align_selected_and_tree.add_argument('--threads',
                                                type=int,
                                                default=1,
                                                help='Number of threads to use for each concurrent alignment. Default '
                                                     'is: %(default)s')
    parser_align_selected_and_tree.add_argument('--no_stitched_contigs',
                                                action='store_true',
                                                default=False,
                                                help='If specified, realign mafft alignments with clustal omega. '
                                                     'Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--use_muscle',
                                                action='store_true',
                                                default=False,
                                                help='If specified, use muscle rather than mafft for initial '
                                                     'alignments. Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--mafft_algorithm',
                                                default='auto',
                                                help='Algorithm to use for mafft alignments. Default is: %(default)s')

    return parser_align_selected_and_tree
