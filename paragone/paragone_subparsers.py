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
                                        help='Do not trim alignments using Trimal. Default is: %(default)s')
    parser_check_and_align.add_argument('--trimal_terminalonly_off',
                                        action='store_true',
                                        default=False,
                                        help='Consider all alignment positions when trimming using Trimal, '
                                             'rather than only terminal columns. Default is: %(default)s')
    parser_check_and_align.add_argument('--trimal_gapthreshold',
                                        type=float,
                                        default=0.12,
                                        help='1 - (fraction of sequences with a gap allowed) when trimming alignments '
                                             'with Trimal. Range: [0 - 1]. Default is: %(default)s')
    parser_check_and_align.add_argument('--trimal_simthreshold',
                                        type=float,
                                        help='Trimal minimum average similarity allowed when trimming alignments '
                                             'with Trimal. Range: [0 - 1]')
    parser_check_and_align.add_argument('--trimal_cons',
                                        type=int,
                                        help='Minimum percentage of positions in the original alignment to conserve '
                                             'when trimming alignments with Trimal. Range: [0 - 100].')
    parser_check_and_align.add_argument('--trimal_nogaps',
                                        action='store_true',
                                        default=False,
                                        help='Remove all positions with gaps in the alignment when trimming. Default '
                                             'is: %(default)s')
    parser_check_and_align.add_argument('--trimal_noallgaps',
                                        action='store_true',
                                        default=False,
                                        help='Remove columns composed only by gaps when trimming alignments. Default '
                                             'is: %(default)s')
    group_1 = parser_check_and_align.add_mutually_exclusive_group(required=False)
    group_1.add_argument('--trimal_gappyout',
                         action='store_const',
                         const='gappyout',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "gappyout" mode. '
                              'This method only uses information based on gaps distribution')
    group_1.add_argument('--trimal_strict',
                         action='store_const',
                         const='strict',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "strict" mode.')
    group_1.add_argument('--trimal_strictplus',
                         action='store_const',
                         const='strictplus',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "strictplus" mode.')
    group_1.add_argument('--trimal_automated1',
                         action='store_const',
                         const='automated1',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal,use a heuristic selection of the automatic method '
                              'based on similarity statistics.')
    parser_check_and_align.add_argument('--trimal_block',
                                        type=int,
                                        help='Minimum column block size to be kept in the trimmed alignment. Available '
                                             'with manual and automatic (gappyout) methods when trimming alignments '
                                             'with Trimal.')
    parser_check_and_align.add_argument('--trimal_resoverlap',
                                        type=float,
                                        help='Minimum overlap of a positions with other positions in the column to be '
                                             'considered a "good position" when trimming alignments with Trimal. '
                                             'Range: [0 - 1].')
    parser_check_and_align.add_argument('--trimal_seqoverlap',
                                        type=int,
                                        help='Minimum percentage of "good positions" that a sequence must have in '
                                             'order to be conserved when trimming alignments with Trimal. '
                                             'Range: [0 - 100].')
    parser_check_and_align.add_argument('--trimal_w',
                                        type=int,
                                        help='(half) Window size, score of position i is the average of the window '
                                             '(i - n) to (i + n), when trimming alignments with Trimal.')
    parser_check_and_align.add_argument('--trimal_gw',
                                        type=int,
                                        default=1,
                                        help='(half) Window size only applies to statistics/methods based on gaps, '
                                             'when trimming alignments with Trimal.')
    parser_check_and_align.add_argument('--trimal_sw',
                                        type=int,
                                        help='(half) Window size only applies to statistics/methods based on '
                                             'similarity, when trimming alignments with Trimal.')
    parser_check_and_align.add_argument('--no_cleaning',
                                        action='store_true',
                                        default=False,
                                        help='Do not clean alignments using TAPER. Default is: %(default)s')
    parser_check_and_align.add_argument('--cleaning_cutoff',
                                        type=int,
                                        default=3,
                                        help='Cutoff value to pass to TAPER. Lower will perform more aggressive '
                                             'cleaning. Default is: %(default)s')
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
    parser_qc_trees_and_fasta.add_argument('mask_tips_alignment_directory',
                                           type=str,
                                           help='directory containing original fasta alignment files')
    parser_qc_trees_and_fasta.add_argument('--mask_tips_alignment_file_suffix',
                                           type=str,
                                           default='.fasta',
                                           help='Suffix for alignment files. Default is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('--min_tips',
                                           type=int,
                                           default=4,
                                           help='The minimum number of tips in a tree after trimming/masking tips or '
                                                'pruning deep paralogs; if below this value, no output tree is '
                                                'written. Default is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('--treeshrink_q_value',
                                           type=float,
                                           default=0.05,
                                           help='q value for TreeShrink; the quantile(s) to set threshold. Default '
                                                'is: %(default)s')
    parser_qc_trees_and_fasta.add_argument('--cut_deep_paralogs_internal_branch_length_cutoff',
                                           type=float,
                                           default=0.3,
                                           help='Internal branch length cutoff cutting tree. Default is: %(default)s')
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
                                                help='Do not trim alignments using Trimal. Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--trimal_terminalonly_off',
                                                action='store_true',
                                                default=False,
                                                help='Consider all alignment positions when trimming using Trimal, '
                                                     'rather than only terminal columns. Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--trimal_gapthreshold',
                                                type=float,
                                                default=0.12,
                                                help='1 - (fraction of sequences with a gap allowed) when trimming '
                                                     'alignments with Trimal. Range: [0 - 1]. Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--trimal_simthreshold',
                                                type=float,
                                                help='Trimal minimum average similarity allowed when trimming '
                                                     'alignments with Trimal. Range: [0 - 1]')
    parser_align_selected_and_tree.add_argument('--trimal_cons',
                                                type=int,
                                                help='Minimum percentage of positions in the original alignment to '
                                                     'conserve when trimming alignments with Trimal. Range: [0 - 100].')
    parser_align_selected_and_tree.add_argument('--trimal_nogaps',
                                                action='store_true',
                                                default=False,
                                                help='Remove all positions with gaps in the alignment when trimming. '
                                                     'Default is: %(default)s')
    parser_align_selected_and_tree.add_argument('--trimal_noallgaps',
                                                action='store_true',
                                                default=False,
                                                help='Remove columns composed only by gaps when trimming alignments. '
                                                     'Default is: %(default)s')
    group_1 = parser_align_selected_and_tree.add_mutually_exclusive_group(required=False)
    group_1.add_argument('--trimal_gappyout',
                         action='store_const',
                         const='gappyout',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "gappyout" mode. '
                              'This method only uses information based on gaps distribution')
    group_1.add_argument('--trimal_strict',
                         action='store_const',
                         const='strict',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "strict" mode.')
    group_1.add_argument('--trimal_strictplus',
                         action='store_const',
                         const='strictplus',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "strictplus" mode.')
    group_1.add_argument('--trimal_automated1',
                         action='store_const',
                         const='automated1',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal,use a heuristic selection of the automatic method '
                              'based on similarity statistics.')
    parser_align_selected_and_tree.add_argument('--trimal_block',
                                                type=int,
                                                help='Minimum column block size to be kept in the trimmed alignment. '
                                                     'Available with manual and automatic (gappyout) methods when '
                                                     'trimming alignments with Trimal.')
    parser_align_selected_and_tree.add_argument('--trimal_resoverlap',
                                                type=float,
                                                help='Minimum overlap of a positions with other positions in the '
                                                     'column to be considered a "good position" when trimming '
                                                     'alignments with Trimal. Range: [0 - 1].')
    parser_align_selected_and_tree.add_argument('--trimal_seqoverlap',
                                                type=int,
                                                help='Minimum percentage of "good positions" that a sequence must '
                                                     'have in order to be conserved when trimming alignments with '
                                                     'Trimal. Range: [0 - 100].')
    parser_align_selected_and_tree.add_argument('--trimal_w',
                                                type=int,
                                                help='(half) Window size, score of position i is the average of the '
                                                     'window (i - n) to (i + n), when trimming alignments with Trimal.')
    parser_align_selected_and_tree.add_argument('--trimal_gw',
                                                type=int,
                                                default=1,
                                                help='(half) Window size only applies to statistics/methods based on '
                                                     'gaps, when trimming alignments with Trimal.')
    parser_align_selected_and_tree.add_argument('--trimal_sw',
                                                type=int,
                                                help='(half) Window size only applies to statistics/methods based on '
                                                     'similarity, when trimming alignments with Trimal.')
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
    parser_prune_paralogs.add_argument('--mo_algorithm_paragone',
                                       action='store_true',
                                       default=False,
                                       help='If pruning trees using the MO algorithm, use an updated ParaGone '
                                            'implementation rather than the original Yang and '
                                            'Smith 2014 implementation. Default is: %(default)s')
    parser_prune_paralogs.add_argument('--minimum_taxa',
                                       type=int,
                                       default=4,
                                       help='Minimum number of taxa required. Default is: %(default)s')
    parser_prune_paralogs.add_argument('--ignore_1to1_orthologs',
                                       action='store_true',
                                       default=False,
                                       help='Do not output 1to1 orthologs, i.e. trees with no paralogs. Default is: %('
                                            'default)s')
    parser_prune_paralogs.add_argument('--run_profiler',
                                       action='store_true',
                                       dest='run_profiler',
                                       default=False,
                                       help='If supplied, run the subcommand using cProfile. Saves a *.csv file of '
                                            'results')
    parser_prune_paralogs.add_argument('--debug',
                                       action='store_true',
                                       default=False,
                                       help='If supplied, log additional information when running the subcommand. '
                                            'This can make the log files much larger.')

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
                                         help='Do not trim alignments using Trimal. Default is: %(default)s')
    parser_final_alignments.add_argument('--trimal_terminalonly_off',
                                         action='store_true',
                                         default=False,
                                         help='Consider all alignment positions when trimming using Trimal, '
                                              'rather than only terminal columns. Default is: %(default)s')
    parser_final_alignments.add_argument('--trimal_gapthreshold',
                                         type=float,
                                         default=0.12,
                                         help='1 - (fraction of sequences with a gap allowed) when trimming alignments '
                                              'with Trimal. Range: [0 - 1]. Default is: %(default)s')
    parser_final_alignments.add_argument('--trimal_simthreshold',
                                         type=float,
                                         help='Trimal minimum average similarity allowed when trimming alignments '
                                              'with Trimal. Range: [0 - 1]')
    parser_final_alignments.add_argument('--trimal_cons',
                                         type=int,
                                         help='Minimum percentage of positions in the original alignment to conserve '
                                              'when trimming alignments with Trimal. Range: [0 - 100].')
    parser_final_alignments.add_argument('--trimal_nogaps',
                                         action='store_true',
                                         default=False,
                                         help='Remove all positions with gaps in the alignment when trimming. Default '
                                              'is: %(default)s')
    parser_final_alignments.add_argument('--trimal_noallgaps',
                                         action='store_true',
                                         default=False,
                                         help='Remove columns composed only by gaps when trimming alignments. Default '
                                              'is: %(default)s')
    group_1 = parser_final_alignments.add_mutually_exclusive_group(required=False)
    group_1.add_argument('--trimal_gappyout',
                         action='store_const',
                         const='gappyout',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "gappyout" mode. '
                              'This method only uses information based on gaps distribution')
    group_1.add_argument('--trimal_strict',
                         action='store_const',
                         const='strict',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "strict" mode.')
    group_1.add_argument('--trimal_strictplus',
                         action='store_const',
                         const='strictplus',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "strictplus" mode.')
    group_1.add_argument('--trimal_automated1',
                         action='store_const',
                         const='automated1',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal,use a heuristic selection of the automatic method '
                              'based on similarity statistics.')
    parser_final_alignments.add_argument('--trimal_block',
                                         type=int,
                                         help='Minimum column block size to be kept in the trimmed alignment. '
                                              'Available with manual and automatic (gappyout) methods when trimming '
                                              'alignments with Trimal.')
    parser_final_alignments.add_argument('--trimal_resoverlap',
                                         type=float,
                                         help='Minimum overlap of a positions with other positions in the column to be '
                                              'considered a "good position" when trimming alignments with Trimal. '
                                              'Range: [0 - 1].')
    parser_final_alignments.add_argument('--trimal_seqoverlap',
                                         type=int,
                                         help='Minimum percentage of "good positions" that a sequence must have in '
                                              'order to be conserved when trimming alignments with Trimal. '
                                              'Range: [0 - 100].')
    parser_final_alignments.add_argument('--trimal_w',
                                         type=int,
                                         help='(half) Window size, score of position i is the average of the window '
                                              '(i - n) to (i + n), when trimming alignments with Trimal.')
    parser_final_alignments.add_argument('--trimal_gw',
                                         type=int,
                                         default=1,
                                         help='(half) Window size only applies to statistics/methods based on gaps, '
                                              'when trimming alignments with Trimal.')
    parser_final_alignments.add_argument('--trimal_sw',
                                         type=int,
                                         help='(half) Window size only applies to statistics/methods based on '
                                              'similarity, when trimming alignments with Trimal.')
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
                                      help='Do not trim alignments using Trimal. Default is: %(default)s')
    parser_full_pipeline.add_argument('--trimal_terminalonly_off',
                                      action='store_true',
                                      default=False,
                                      help='Consider all alignment positions when trimming using Trimal, '
                                           'rather than only terminal columns. Default is: %(default)s')
    parser_full_pipeline.add_argument('--trimal_gapthreshold',
                                      type=float,
                                      default=0.12,
                                      help='1 - (fraction of sequences with a gap allowed) when trimming alignments '
                                           'with Trimal. Range: [0 - 1]. Default is: %(default)s')
    parser_full_pipeline.add_argument('--trimal_simthreshold',
                                      type=float,
                                      help='Trimal minimum average similarity allowed when trimming alignments '
                                           'with Trimal. Range: [0 - 1]')
    parser_full_pipeline.add_argument('--trimal_cons',
                                      type=int,
                                      help='Minimum percentage of positions in the original alignment to conserve '
                                           'when trimming alignments with Trimal. Range: [0 - 100].')
    parser_full_pipeline.add_argument('--trimal_nogaps',
                                      action='store_true',
                                      default=False,
                                      help='Remove all positions with gaps in the alignment when trimming. Default '
                                           'is: %(default)s')
    parser_full_pipeline.add_argument('--trimal_noallgaps',
                                      action='store_true',
                                      default=False,
                                      help='Remove columns composed only by gaps when trimming alignments. Default '
                                           'is: %(default)s')
    group_1 = parser_full_pipeline.add_mutually_exclusive_group(required=False)
    group_1.add_argument('--trimal_gappyout',
                         action='store_const',
                         const='gappyout',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "gappyout" mode. '
                              'This method only uses information based on gaps distribution')
    group_1.add_argument('--trimal_strict',
                         action='store_const',
                         const='strict',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "strict" mode.')
    group_1.add_argument('--trimal_strictplus',
                         action='store_const',
                         const='strictplus',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal, use the automated selection on "strictplus" mode.')
    group_1.add_argument('--trimal_automated1',
                         action='store_const',
                         const='automated1',
                         dest='automated_method',
                         default=False,
                         help='When trimming alignments with Trimal,use a heuristic selection of the automatic method '
                              'based on similarity statistics.')
    parser_full_pipeline.add_argument('--trimal_block',
                                      type=int,
                                      help='Minimum column block size to be kept in the trimmed alignment. '
                                           'Available with manual and automatic (gappyout) methods when trimming '
                                           'alignments with Trimal.')
    parser_full_pipeline.add_argument('--trimal_resoverlap',
                                      type=float,
                                      help='Minimum overlap of a positions with other positions in the column to be '
                                           'considered a "good position" when trimming alignments with Trimal. '
                                           'Range: [0 - 1].')
    parser_full_pipeline.add_argument('--trimal_seqoverlap',
                                      type=int,
                                      help='Minimum percentage of "good positions" that a sequence must have in order '
                                           'to be conserved when trimming alignments with Trimal. Range: [0 - 100].')
    parser_full_pipeline.add_argument('--trimal_w',
                                      type=int,
                                      help='(half) Window size, score of position i is the average of the window '
                                           '(i - n) to (i + n), when trimming alignments with Trimal.')
    parser_full_pipeline.add_argument('--trimal_gw',
                                      type=int,
                                      default=1,
                                      help='(half) Window size only applies to statistics/methods based on gaps, '
                                           'when trimming alignments with Trimal.')
    parser_full_pipeline.add_argument('--trimal_sw',
                                      type=int,
                                      help='(half) Window size only applies to statistics/methods based on '
                                           'similarity, when trimming alignments with Trimal.')
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
                                      help='Do not clean alignments using TAPER. Default is: %(default)s')
    parser_full_pipeline.add_argument('--cleaning_cutoff',
                                      type=int,
                                      default=3,
                                      help='Cutoff value to pass to TAPER. Lower will perform more aggressive '
                                           'cleaning. Default is: %(default)s')
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
                                      default=0.05,
                                      help='q value for TreeShrink; the quantile(s) to set threshold. Default '
                                           'is: %(default)s')
    parser_full_pipeline.add_argument('--cut_deep_paralogs_internal_branch_length_cutoff',
                                      type=float,
                                      default=0.3,
                                      help='Internal branch length cutoff cutting tree. Default is: '
                                           '%(default)s')
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
    parser_full_pipeline.add_argument('--mo_algorithm_paragone',
                                      action='store_true',
                                      default=False,
                                      help='If pruning trees using the MO algorithm, use an updated ParaGone '
                                           'implementation rather than the original Yang and '
                                           'Smith 2014 implementation. Default is: %(default)s')
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
    parser_full_pipeline.add_argument('--debug',
                                      action='store_true',
                                      default=False,
                                      help='If supplied, log additional information when running the subcommand. '
                                           'This can make the log files much larger.')
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
