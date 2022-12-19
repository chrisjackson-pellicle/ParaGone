#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au

"""
Contains argument subparsers
"""


def add_check_and_align_parser(subparsers):
    """
    Parser for check_and_align

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_check_and_align = subparsers.add_parser('check_and_align',
                                                   help='Check input files and outgroup coverage; generate per-gene '
                                                        'paralog alignments; clean alignments')
    parser_check_and_align.add_argument('gene_fasta_directory',
                                        type=str,
                                        help='Directory contains fasta files with paralog sequences')
    parser_check_and_align.add_argument('--gene_name_delimiter',
                                        type=str,
                                        default='_',
                                        help='Delimiter in paralog filename to extract gene name. Default is: '
                                             '%(default)s')
    parser_check_and_align.add_argument('--gene_name_field_num',
                                        type=int,
                                        default='1',
                                        help='From paralog filename, number of fields to extract gene name. Default '
                                             'is: %(default)s')
    parser_check_and_align.add_argument('--external_outgroups_file',
                                        type=str,
                                        default=None,
                                        help='File in fasta format with additional outgroup sequences to add to each '
                                             'gene')
    parser_check_and_align.add_argument('--external_outgroup',
                                        action='append',
                                        type=str,
                                        dest='external_outgroups',
                                        default=None,
                                        help='If a taxon name is provided, only use these sequences from '
                                             'the user-provided external_outgroups_file. Note that this parameter can '
                                             'be specified one ore more times.')
    parser_check_and_align.add_argument('--internal_outgroup',
                                        action='append',
                                        type=str,
                                        dest='internal_outgroups',
                                        default=None,
                                        help='Taxon name to use as an internal outgroup (i.e. present in input '
                                             'paralog files). Note that this parameter can be specified one or more '
                                             'times.')
    parser_check_and_align.add_argument('--pool',
                                        type=int,
                                        default=1,
                                        help='Number of alignments to run concurrently. Default is: %(default)s')
    parser_check_and_align.add_argument('--threads',
                                        type=int,
                                        default=1,
                                        help='Number of threads to use for each concurrent alignment. Default '
                                             'is: %(default)s')
    parser_check_and_align.add_argument('--no_stitched_contigs',
                                        action='store_true',
                                        default=False,
                                        help='If specified, realign mafft alignments with clustal omega. Default is: '
                                             '%(default)s')
    parser_check_and_align.add_argument('--use_muscle',
                                        action='store_true',
                                        default=False,
                                        help='If specified, use muscle rather than mafft for initial alignments. '
                                             'Default is: %(default)s')
    parser_check_and_align.add_argument('--mafft_algorithm',
                                        default='auto',
                                        help='Algorithm to use for mafft alignments. Default is: %(default)s')
    parser_check_and_align.add_argument('--no_trimming',
                                        action='store_true',
                                        default=False,
                                        help='No not trim alignments using Trimal. Default is: %(default)s')
    parser_check_and_align.add_argument('--no_cleaning',
                                        action='store_true',
                                        default=False,
                                        help='No not clean alignments using HmmCleaner.pl. Default is: %(default)s')
    parser_check_and_align.add_argument('--run_profiler',
                                        action='store_true',
                                        dest='run_profiler',
                                        default=False,
                                        help='If supplied, run the subcommand using cProfile. Saves a *.csv file of '
                                             'results')

    return parser_check_and_align


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
    parser_alignment_to_tree.add_argument('--run_profiler',
                                          action='store_true',
                                          dest='run_profiler',
                                          default=False,
                                          help='If supplied, run the subcommand using cProfile. Saves a *.csv file of '
                                               'results')

    return parser_alignment_to_tree

#
# def add_collate_alignments_and_trees_parser(subparsers):
#     """
#     Parser for collate_alignments_and_trees
#
#     :param argparse._SubParsersAction subparsers:
#     :return None: no return value specified; default is None
#     """
#
#     parser_collate_alignments_and_trees = subparsers.add_parser(
#         'collate_alignments_and_trees',
#         help='Collates all HmmCleaned alignments into a single folder. Collates all corresponding trees into a '
#              'single folder')
#     group_1 = parser_collate_alignments_and_trees.add_mutually_exclusive_group(required=True)
#     group_1.add_argument('--from_alignment_to_tree',
#                          action='store_true',
#                          dest='from_alignment_to_tree',
#                          default=False,
#                          help='If set, trees are from step "alignment_to_tree".')
#     group_1.add_argument('--from_align_selected_and_tree',
#                          action='store_true',
#                          dest='from_align_selected_and_tree',
#                          default=False,
#                          help='If set, trees are from step "align_selected_and_tree".')
#     group_1.add_argument('--from_prune_paralogs_mo',
#                          action='store_const', const='mo',
#                          dest='from_prune_paralogs',
#                          default=False,
#                          help='If set, sequences are from paralog pruning step "prune_paralogs_mo".')
#     group_1.add_argument('--from_prune_paralogs_rt',
#                          action='store_const', const='rt',
#                          dest='from_prune_paralogs',
#                          default=False,
#                          help='If set, sequences are from paralog pruning step "prune_paralogs_rt".')
#     group_1.add_argument('--from_prune_paralogs_mi',
#                          action='store_const', const='mi',
#                          dest='from_prune_paralogs',
#                          default=False,
#                          help='If set, sequences are from paralog pruning step "prune_paralogs_mi".')
#
#     parser_collate_alignments_and_trees.add_argument('--tree_file_suffix',
#                                                      type=str,
#                                                      default='.treefile',
#                                                      help='Suffix for newick tree files. Default is: %(default)s')
#     parser_collate_alignments_and_trees.add_argument('--run_profiler',
#                                                      action='store_true',
#                                                      dest='run_profiler',
#                                                      default=False,
#                                                      help='If supplied, run the subcommand using cProfile. Saves a '
#                                                           '*.csv file of results')
#
#     return parser_collate_alignments_and_trees


def add_qc_trees_and_extract_fasta(subparsers):
    """
    Parser for add_qc_trees_and_fasta

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_qc_trees_and_fasta = subparsers.add_parser('qc_trees_and_extract_fasta',
                                                      help='Quality control trees; for remaining tips, '
                                                           'extract corresponding fasta sequences')
    # parser_qc_trees_and_fasta.add_argument('treefile_directory',
    #                                        type=str,
    #                                        help='directory containing tree newick files')
    # parser_qc_trees_and_fasta.add_argument('--tree_file_suffix',
    #                                        type=str,
    #                                        default='.treefile',
    #                                        help='Suffix for newick tree files. Default is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('--trim_tips_relative_cutoff',
                                           type=float,
                                           default=0.2,
                                           help='Relative cutoff for removing tree tips. Default is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('--trim_tips_absolute_cutoff',
                                           type=float,
                                           default=0.4,
                                           help='Absolute cutoff for removing tree tips. Default is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('mask_tips_alignment_directory',
                                           type=str,
                                           help='directory containing original fasta alignment files')
    parser_qc_trees_and_fasta.add_argument('--mask_tips_alignment_file_suffix',
                                           type=str,
                                           default='.fasta',
                                           help='Suffix for alignment files. Default is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('--mask_tips_remove_paraphyletic_tips',
                                           action='store_true',
                                           default=False,
                                           help='Remove paraphyletic tree tips. Default is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('--cut_deep_paralogs_internal_branch_length_cutoff',
                                           type=float,
                                           default=0.3,
                                           help='Internal branch length cutoff cutting tree. Default is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('--cut_deep_paralogs_minimum_number_taxa',
                                           type=int,
                                           default=4,
                                           help='Minimum number of taxa in tree for tree to be retained. Default is: '
                                                '%(default)s')
    parser_qc_trees_and_fasta.add_argument('--run_profiler',
                                           action='store_true',
                                           dest='run_profiler',
                                           default=False,
                                           help='If supplied, run the subcommand using cProfile. Saves a '
                                                '*.csv file of results')

    # Set defaults for subparser <parser_qc_trees_and_fasta>:
    parser_qc_trees_and_fasta.set_defaults(
        treefile_directory='07_trees_pre_quality_control_trimmed_masked_cut',
        # alignment_directory='03_input_paralog_fasta_with_sanitised_filenames_alignments_hmmcleaned',
        # tree_file_suffix='.subtree',
        from_cut_deep_paralogs=True)

    return parser_qc_trees_and_fasta


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
    group_1 = parser_fasta_from_tree.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--from_cut_deep_paralogs',
                         action='store_true',
                         dest='from_cut_deep_paralogs',
                         default=False,
                         help='If set, trees are from QC step "cut_deep_paralogs".')
    group_1.add_argument('--from_prune_paralogs_mo',
                         action='store_const', const='mo',
                         dest='from_prune_paralogs',
                         default=False,
                         help='If set, trees are from paralog pruning step "prune_paralogs_mo".')
    group_1.add_argument('--from_prune_paralogs_rt',
                         action='store_const', const='rt',
                         dest='from_prune_paralogs',
                         default=False,
                         help='If set, trees are from paralog pruning step ="prune_paralogs_rt".')
    group_1.add_argument('--from_prune_paralogs_mi',
                         action='store_const', const='mi',
                         dest='from_prune_paralogs',
                         default=False,
                         help='If set, trees are from paralog pruning step "prune_paralogs_mi".')
    # parser_fasta_from_tree.add_argument('--batch_size',
    #                                     type=int,
    #                                     default=20,
    #                                     help='Number of fasta files in each batch, from input paralog fasta files. '
    #                                          'Default is: %(default)s')
    # parser_fasta_from_tree.add_argument('--run_profiler',
    #                                     action='store_true',
    #                                     dest='run_profiler',
    #                                     default=False,
    #                                     help='If supplied, run the subcommand using cProfile. Saves a '
    #                                          '*.csv file of results')

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
    parser_align_selected_and_tree.add_argument('--run_profiler',
                                                action='store_true',
                                                dest='run_profiler',
                                                default=False,
                                                help='If supplied, run the subcommand using cProfile. Saves a '
                                                     '*.csv file of results')

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
    parser_prune_paralogs_mo.add_argument('--run_profiler',
                                          action='store_true',
                                          dest='run_profiler',
                                          default=False,
                                          help='If supplied, run the subcommand using cProfile. Saves a *.csv file of '
                                               'results')

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
    parser_prune_paralogs_rt.add_argument('--run_profiler',
                                          action='store_true',
                                          dest='run_profiler',
                                          default=False,
                                          help='If supplied, run the subcommand using cProfile. Saves a *.csv file of '
                                               'results')

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
    parser_prune_paralogs_mi.add_argument('--run_profiler',
                                          action='store_true',
                                          dest='run_profiler',
                                          default=False,
                                          help='If supplied, run the subcommand using cProfile. Saves a *.csv file of '
                                               'results')

    return parser_prune_paralogs_mi


def add_strip_names_and_align_parser(subparsers):
    """
    Parser for strip_names_and_align

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_align_selected_and_tree = subparsers.add_parser('strip_names_and_align',
                                                           help='Strip names of paralog designations (e.g. .main, .0, '
                                                                '.1 etc), and performs a final alignment step')
    parser_align_selected_and_tree.add_argument('selected_alignment_directory',
                                                type=str,
                                                help='directory containing selected alignment files corresponding to '
                                                     'pruned trees from one of MO, RT, MI algorithms.')
    group_1 = parser_align_selected_and_tree.add_mutually_exclusive_group(required=True)
    group_1.add_argument('--from_prune_paralogs_mo',
                         action='store_const', const='mo',
                         dest='from_prune_paralogs',
                         default=False,
                         help='If set, sequences are from paralog pruning step "prune_paralogs_mo".')
    group_1.add_argument('--from_prune_paralogs_rt',
                         action='store_const', const='rt',
                         dest='from_prune_paralogs',
                         default=False,
                         help='If set, sequences are from paralog pruning step ="prune_paralogs_rt".')
    group_1.add_argument('--from_prune_paralogs_mi',
                         action='store_const', const='mi',
                         dest='from_prune_paralogs',
                         default=False,
                         help='If set, sequences are from paralog pruning step "prune_paralogs_mi".')
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
    parser_align_selected_and_tree.add_argument('--run_profiler',
                                                action='store_true',
                                                dest='run_profiler',
                                                default=False,
                                                help='If supplied, run the subcommand using cProfile. Saves a *.csv '
                                                     'file of results')

    return parser_align_selected_and_tree
