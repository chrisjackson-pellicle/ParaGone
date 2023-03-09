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
                                                   help='Check input files and outgroup coverage;\ngenerate per-gene '
                                                        'paralog alignments;\ntrim and/or clean alignments (optional)')
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
    parser_check_and_align.add_argument('--use_clustal',
                                        action='store_true',
                                        default=False,
                                        help='If specified, alignments are performed using Clustal Omega rather than '
                                             'MAFFT. If the flag "--mafft_adjustdirection" is also provided, '
                                             'alignments are performed with MAFFT first, followed by realignment '
                                             'using with Clustal Omega. Default is: %(default)s')
    parser_check_and_align.add_argument('--mafft_algorithm',
                                        default='auto',
                                        help='Algorithm to use for MAFFT alignments. Default is: %(default)s')
    parser_check_and_align.add_argument('--mafft_adjustdirection',
                                        action='store_true',
                                        default=False,
                                        help='Allow MAFFT to generate reverse complement sequences, as necessary, '
                                             'and align them together with the remaining sequences. Note that '
                                             'the first sequence is assumed to be in the correct orientation. Default '
                                             'is: %(default)s')
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
                                                     help='Generate phylogenetic trees from alignments')
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


def add_qc_trees_and_extract_fasta(subparsers):
    """
    Parser for add_qc_trees_and_fasta

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_qc_trees_and_fasta = subparsers.add_parser('qc_trees_and_extract_fasta',
                                                      help='Quality control trees;\nextract fasta sequences '
                                                           'corresponding to remaining tips')
    parser_qc_trees_and_fasta.add_argument('--min_tips',
                                           type=int,
                                           default=4,
                                           help='The minimum number of tips in a tree after trimming/masking tips or '
                                                'pruning deep paralogs; if below this value, no output tree is '
                                                'written. Default is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('--treeshrink_q_value',
                                           type=float,
                                           default=0.20,
                                           help='q value for TreeShrink; the quantile(s) to set threshold. Default '
                                                'is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('mask_tips_alignment_directory',
                                           type=str,
                                           help='directory containing original fasta alignment files')
    parser_qc_trees_and_fasta.add_argument('--mask_tips_alignment_file_suffix',
                                           type=str,
                                           default='.fasta',
                                           help='Suffix for alignment files. Default is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('--cut_deep_paralogs_internal_branch_length_cutoff',
                                           type=float,
                                           default=0.3,
                                           help='Internal branch length cutoff cutting tree. Default is: %(default)s')
    # parser_qc_trees_and_fasta.add_argument('--cut_deep_paralogs_minimum_number_taxa',
    #                                        type=int,
    #                                        default=4,
    #                                        help='Minimum number of taxa in tree for tree to be retained. Default is: '
    #                                             '%(default)s')
    parser_qc_trees_and_fasta.add_argument('--run_profiler',
                                           action='store_true',
                                           dest='run_profiler',
                                           default=False,
                                           help='If supplied, run the subcommand using cProfile. Saves a '
                                                '*.csv file of results')

    # Set defaults for subparser <parser_qc_trees_and_fasta>:
    parser_qc_trees_and_fasta.set_defaults(
        from_cut_deep_paralogs=True)

    return parser_qc_trees_and_fasta


def add_align_selected_and_tree_parser(subparsers):
    """
    Parser for align_selected_and_tree

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_align_selected_and_tree = subparsers.add_parser('align_selected_and_tree',
                                                           help='Align selected fasta seqs for each subtree, '
                                                                'and generate a new tree')
    parser_align_selected_and_tree.add_argument('selected_alignment_directory',
                                                type=str,
                                                help='directory containing selected alignment files corresponding to '
                                                     'subtrees')
    parser_align_selected_and_tree.add_argument('qc_alignment_directory',
                                                type=str,
                                                help='directory containing quality controlled alignment files')
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
    parser_align_selected_and_tree.add_argument('--use_clustal',
                                                action='store_true',
                                                default=False,
                                                help='If specified, alignments are performed using Clustal Omega '
                                                     'rather than MAFFT. If the flag "--mafft_adjustdirection" is '
                                                     'also provided, alignments are performed with MAFFT first, '
                                                     'followed by realignment using with Clustal Omega. Default is: '
                                                     '%(default)s')
    parser_align_selected_and_tree.add_argument('--mafft_algorithm',
                                                default='auto',
                                                help='Algorithm to use for MAFFT alignments. Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--no_trimming',
                                                action='store_true',
                                                default=False,
                                                help='No not trim alignments using Trimal. Default is: %(default)s')
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


def add_prune_paralogs_parser(subparsers):
    """
    Parser for prune_paralogs

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_prune_paralogs = subparsers.add_parser('prune_paralogs',
                                                  help='Prune paralogs from tree using one or more of the '
                                                       'algorithms: \nMonophyletic Outgroups (MO), RooTed ingroups ('
                                                       'RT), Maximum Inclusion (MI)')
    parser_prune_paralogs.add_argument('--mo',
                                       action='store_true',
                                       default=False,
                                       help='Run the Monophyletic Outgroups (MO) algorithm')
    parser_prune_paralogs.add_argument('--mi',
                                       action='store_true',
                                       default=False,
                                       help='Run the Maximum Inclusion (MI) algorithm')
    parser_prune_paralogs.add_argument('--rt',
                                       action='store_true',
                                       default=False,
                                       help='Run the RooTed ingroups (RT) algorithm')
    parser_prune_paralogs.add_argument('--minimum_taxa',
                                       type=int,
                                       default=4,
                                       help='Minimum number of taxa required. Default is: %(default)s')
    parser_prune_paralogs.add_argument('--ignore_1to1_orthologs',
                                       action='store_true',
                                       default=False,
                                       help='Output 1to1 orthologs, i.e. trees with no paralogs. Default is: %('
                                            'default)s')
    # parser_prune_paralogs.add_argument('--relative_tip_cutoff',
    #                                    type=float,
    #                                    default=0.2,
    #                                    help='Relative tip cut-off threshold. Default is: %(default)s')
    # parser_prune_paralogs.add_argument('--absolute_tip_cutoff',
    #                                    type=float,
    #                                    default=0.4,
    #                                    help='Absolute tip cut-off threshold. Default is: %(default)s')
    parser_prune_paralogs.add_argument('--run_profiler',
                                       action='store_true',
                                       dest='run_profiler',
                                       default=False,
                                       help='If supplied, run the subcommand using cProfile. Saves a *.csv file of '
                                            'results')

    return parser_prune_paralogs


def add_final_alignments_parser(subparsers):
    """
    Parser for final_alignments

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_final_alignments = subparsers.add_parser('final_alignments',
                                                    help='Recover fasta sequences for pruned trees;\nstrip names of '
                                                         'paralog designations (e.g. .main, .0, .1 etc);\nperform '
                                                         'final alignments')
    parser_final_alignments.add_argument('--mo',
                                         action='store_true',
                                         default=False,
                                         help='Run the Monophyletic Outgroups (MO) algorithm')
    parser_final_alignments.add_argument('--mi',
                                         action='store_true',
                                         default=False,
                                         help='Run the Maximum Inclusion (MI) algorithm')
    parser_final_alignments.add_argument('--rt',
                                         action='store_true',
                                         default=False,
                                         help='Run the RooTed ingroups (RT) algorithm')
    parser_final_alignments.add_argument('--pool',
                                         type=int,
                                         default=1,
                                         help='Number of alignments to run concurrently. Default is: %('
                                              'default)s')
    parser_final_alignments.add_argument('--threads',
                                         type=int,
                                         default=1,
                                         help='Number of threads to use for each concurrent alignment. Default '
                                              'is: %(default)s')
    parser_final_alignments.add_argument('--use_clustal',
                                         action='store_true',
                                         default=False,
                                         help='If specified, alignments are performed using Clustal Omega rather than '
                                              'MAFFT. If the flag "--mafft_adjustdirection" is also provided, '
                                              'alignments are performed with MAFFT first, followed by realignment '
                                              'using with Clustal Omega. Default is: %(default)s')
    parser_final_alignments.add_argument('--mafft_algorithm',
                                         default='auto',
                                         help='Algorithm to use for MAFFT alignments. Default is: %(default)s')
    parser_final_alignments.add_argument('--no_trimming',
                                         action='store_true',
                                         default=False,
                                         help='No not trim alignments using Trimal. Default is: %(default)s')
    parser_final_alignments.add_argument('--run_profiler',
                                         action='store_true',
                                         dest='run_profiler',
                                         default=False,
                                         help='If supplied, run the subcommand using cProfile. Saves a *.csv '
                                              'file of results')
    parser_final_alignments.add_argument('--keep_intermediate_files',
                                         action='store_true',
                                         default=False,
                                         help='Keep all intermediate files and folders. Default is: %(default)s')

    parser_final_alignments.set_defaults(
        from_cut_deep_paralogs=False)

    return parser_final_alignments


def add_full_pipeline_parser(subparsers):
    """
    Parser for full_pipeline

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_full_pipeline = subparsers.add_parser('full_pipeline',
                                                 help='Run all steps of the ParaGone pipeline.')
    parser_full_pipeline.add_argument('gene_fasta_directory',
                                      type=str,
                                      help='Directory contains fasta files with paralog sequences')
    parser_full_pipeline.add_argument('--gene_name_delimiter',
                                      type=str,
                                      default='_',
                                      help='Delimiter in paralog filename to extract gene name. Default is: '
                                           '%(default)s')
    parser_full_pipeline.add_argument('--gene_name_field_num',
                                      type=int,
                                      default='1',
                                      help='From paralog filename, number of fields to extract gene name. '
                                           'Default is: %(default)s')
    parser_full_pipeline.add_argument('--external_outgroups_file',
                                      type=str,
                                      default=None,
                                      help='File in fasta format with additional outgroup sequences to add '
                                           'to each gene')
    parser_full_pipeline.add_argument('--external_outgroup',
                                      action='append',
                                      type=str,
                                      dest='external_outgroups',
                                      default=None,
                                      help='If a taxon name is provided, only use these sequences from the '
                                           'user-provided external_outgroups_file. Note that this parameter '
                                           'can be specified one ore more times.')
    parser_full_pipeline.add_argument('--internal_outgroup',
                                      action='append',
                                      type=str,
                                      dest='internal_outgroups',
                                      default=None,
                                      help='Taxon name to use as an internal outgroup (i.e. present in input '
                                           'paralog files). Note that this parameter can be specified one or '
                                           'more times.')
    parser_full_pipeline.add_argument('--pool',
                                      type=int,
                                      default=1,
                                      help='Number of alignments to run concurrently. Default is: %(default)s')
    parser_full_pipeline.add_argument('--threads',
                                      type=int,
                                      default=1,
                                      help='Number of threads to use for each concurrent alignment. Default '
                                           'is: %(default)s')
    parser_full_pipeline.add_argument('--use_clustal',
                                      action='store_true',
                                      default=False,
                                      help='If specified, alignments are performed using Clustal Omega rather than '
                                           'MAFFT. If the flag "--mafft_adjustdirection" is also provided, '
                                           'alignments are performed with MAFFT first, followed by realignment '
                                           'using with Clustal Omega. Default is: %(default)s')
    parser_full_pipeline.add_argument('--mafft_algorithm',
                                      default='auto',
                                      help='Algorithm to use for MAFFT alignments. Default is: %(default)s')
    parser_full_pipeline.add_argument('--no_trimming',
                                      action='store_true',
                                      default=False,
                                      help='No not trim alignments using Trimal. Default is: %(default)s')
    parser_full_pipeline.add_argument('--mafft_adjustdirection',
                                      action='store_true',
                                      default=False,
                                      help='Allow MAFFT to generate reverse complement sequences, as necessary, '
                                           'and align them together with the remaining sequences. Note that '
                                           'the first sequence is assumed to be in the correct orientation. Default '
                                           'is: %(default)s')
    parser_full_pipeline.add_argument('--no_cleaning',
                                      action='store_true',
                                      default=False,
                                      help='No not clean alignments using HmmCleaner.pl. Default is: '
                                           '%(default)s')
    parser_full_pipeline.add_argument('--generate_bootstraps',
                                      action='store_true',
                                      default=False,
                                      help='Create bootstraps for trees using UFBoot. Default is: '
                                           '%(default)s')
    parser_full_pipeline.add_argument('--use_fasttree',
                                      action='store_true',
                                      default=False,
                                      help='Use FastTree instead of IQTREE. Default is: %(default)s')
    parser_full_pipeline.add_argument('--run_profiler',
                                      action='store_true',
                                      dest='run_profiler',
                                      default=False,
                                      help='If supplied, run the subcommand using cProfile. Saves a *.csv '
                                           'file of results')
    parser_full_pipeline.add_argument('--min_tips',
                                      type=int,
                                      default=4,
                                      help='The minimum number of tips in a tree after trimming/masking tips or '
                                           'pruning deep paralogs; if below this value, no output tree is written. '
                                           'Default is: %(default)s')
    parser_full_pipeline.add_argument('--treeshrink_q_value',
                                      type=float,
                                      default=0.20,
                                      help='q value for TreeShrink; the quantile(s) to set threshold. Default '
                                           'is: %(default)s')
    parser_full_pipeline.add_argument('--cut_deep_paralogs_internal_branch_length_cutoff',
                                      type=float,
                                      default=0.3,
                                      help='Internal branch length cutoff cutting tree. Default is: '
                                           '%(default)s')
    # parser_full_pipeline.add_argument('--cut_deep_paralogs_minimum_number_taxa',
    #                                   type=int,
    #                                   default=4,
    #                                   help='Minimum number of taxa in tree for tree to be retained. Default '
    #                                        'is: %(default)s')
    parser_full_pipeline.add_argument('--mo',
                                      action='store_true',
                                      default=False,
                                      help='Run the Monophyletic Outgroups (MO) algorithm')
    parser_full_pipeline.add_argument('--mi',
                                      action='store_true',
                                      default=False,
                                      help='Run the Maximum Inclusion (MI) algorithm')
    parser_full_pipeline.add_argument('--rt',
                                      action='store_true',
                                      default=False,
                                      help='Run the RooTed ingroups (RT) algorithm')
    parser_full_pipeline.add_argument('--minimum_taxa',
                                      type=int,
                                      default=4,
                                      help='Minimum number of taxa required. Default is: %(default)s')
    parser_full_pipeline.add_argument('--ignore_1to1_orthologs',
                                      action='store_true',
                                      default=False,
                                      help='Output 1to1 orthologs, i.e. trees with no paralogs. Default is: '
                                           '%(default)s')
    parser_full_pipeline.add_argument('--keep_intermediate_files',
                                      action='store_true',
                                      default=False,
                                      help='Keep all intermediate files and folders. Default is: %(default)s')
    return parser_full_pipeline


def add_delete_intermediate_files_parser(subparsers):
    """
    Parser for delete_intermediate_files

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_delete_intermediate_files = subparsers.add_parser('delete_intermediate_files',
                                                             help='Deletes all intermediate files and folders after '
                                                                  'the full pipline has been run')
    parser_delete_intermediate_files.add_argument('--run_profiler',
                                                  action='store_true',
                                                  dest='run_profiler',
                                                  default=False,
                                                  help='If supplied, run the subcommand using cProfile. Saves a *.csv '
                                                       'file of results')

    return parser_delete_intermediate_files
