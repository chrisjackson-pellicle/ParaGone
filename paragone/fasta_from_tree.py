# Author: Alexander Schmidt-Lebuhn (modified by Chris Jackson July 2021 chris.jackson@rbg.vic.gov.au)
# https://github.com/chrisjackson-pellicle

"""
Takes a newick tree and a fasta alignment as input. Returns a fasta alignment containing only the sequences
corresponding to the tree tip names.
"""

from Bio import AlignIO
from Bio import Phylo
import Bio.Align
import glob
import os
import sys
import textwrap
import shutil
import re

from paragone import utils


def subsample_alignments(treefile_directory,
                         output_folder,
                         tree_suffix,
                         alignment_directory,
                         trimmed=False,
                         from_cut_deep_paralogs=False,
                         logger=None):
    """
    Takes a pruned/QC'd tree file, finds the original matching alignment, and sub-samples that alignment to recover
    only sequences corresponding to tree tip names.

    :param str treefile_directory: path to directory containing tree newick files
    :param str output_folder: path to output folder for selected fasta alignments
    :param str tree_suffix: suffix for the tree files
    :param str alignment_directory: path to the directory containing fasta alignments
    :param bool trimmed: if True, final pre-resolution tree alignments were trimmed with TrimAl
    :param bool from_cut_deep_paralogs: if True, process tree file names accordingly to recover gene names
    :param logging.Logger logger: a logger object
    :return str, dict output_folder, alignment_filtering_dict: path the output folder with filtered alignments,
    dictionary of filtering stats for each tree/alignment
    """

    logger.info(f'{"[INFO]:":10} Recovering alignment sequences corresponding to tree tip names...')

    # Check for trim/clean status of original alignments:
    if from_cut_deep_paralogs:
        if re.search('cleaned', alignment_directory) and re.search('trimmed', alignment_directory):
            logger.debug(f'Input alignment folder is {alignment_directory}; sequenced were trimmed and cleaned')
            alignment_suffix = '.aln.trimmed.cleaned.fasta'
        elif re.search('cleaned', alignment_directory):
            logger.debug(f'Input alignment folder is {alignment_directory}; sequenced were cleaned but not trimmed')
            alignment_suffix = '.aln.cleaned.fasta'
        elif re.search('trimmed', alignment_directory):
            logger.debug(f'Input alignment folder is {alignment_directory}; sequenced were trimmed but not cleaned')
            alignment_suffix = '.aln.trimmed.fasta'
        else:
            logger.debug(f'Input alignment folder is {alignment_directory}; sequenced were neither trimmed or cleaned')
            alignment_suffix = '.aln.fasta'
    else:
        alignment_suffix = '.outgroup_added.aln.fasta' if not trimmed else '.outgroup_added.aln.trimmed.fasta'

    # Capture number of sequences pre and post filtering in a dictionary for report:
    alignment_filtering_dict = {}

    for tree in glob.glob(f'{treefile_directory}/*{tree_suffix}'):
        read_tree = Phylo.read(tree, "newick")
        tree_terminals = read_tree.get_terminals()
        tree_basename = os.path.basename(tree)

        # Derive the matching alignment file name depending on input tree file name:
        if from_cut_deep_paralogs:  # e.g. 4471_1.subtree
            alignment_prefix = '_'.join(tree_basename.split('_')[0:-1])
            output_alignment_prefix = tree_basename.split('.')[0]
            matching_alignment = f'{alignment_directory}/{alignment_prefix}{alignment_suffix}'
        else:  # e.g. 4691_1.1to1ortho.tre, 4471_1.inclade1.ortho1.tre, 4527_1.MIortho1.tre. etc
            alignment_prefix = tree_basename.split('.')[0]
            output_alignment_prefix = '.'.join(tree_basename.split('.')[0:-1])
            matching_alignment = f'{alignment_directory}/{alignment_prefix}{alignment_suffix}'

        # Read in original alignments and select seqs matching tree termini:
        alignment = AlignIO.read(matching_alignment, "fasta")
        subalignment = Bio.Align.MultipleSeqAlignment([])
        for k in range(0, len(alignment)):
            for j in range(0, len(tree_terminals)):
                if tree_terminals[j].name == alignment[k].id:
                    subalignment.append(alignment[k])

        assert len(tree_terminals) == len(subalignment)

        # Capture data:
        alignment_filtering_dict[tree_basename] = [len(tree_terminals), len(alignment), len(subalignment)]

        # Write an alignment of the sub-selected sequences:
        AlignIO.write(subalignment, f'{output_folder}/{output_alignment_prefix}.selected.fasta', "fasta")

    return output_folder, alignment_filtering_dict


def write_fasta_from_tree_report(alignment_filtering_dict,
                                 report_directory,
                                 from_cut_deep_paralogs,
                                 algorithm_suffix,
                                 logger=None):
    """
    Writes a *.tsv report detailing number of tips in QC'd tree, and number of sequences in original QC'd alignment.

    :param dict alignment_filtering_dict: dictionary of filtering stats for each tree/alignment
    :param str report_directory: path to directory for report files
    :param bool from_cut_deep_paralogs: if True, add 'cut' to report filename
    :param strNone algorithm_suffix: if extracting seqs from pruned trees, the algorithm suffix mo/rt/mi, else None
    :param logging.Logger logger: a logger object
    :return:
    """

    if from_cut_deep_paralogs:
        report_filename = f'{report_directory}/fasta_from_qc_trees_report.tsv'
    elif algorithm_suffix in ['mo', 'mi', 'rt']:
        report_filename = f'{report_directory}/fasta_from_tree_{algorithm_suffix}_report.tsv'

    fill = utils.fill_forward_slash(f'{"[INFO]:":10} Writing fasta_from_tree report to file: "{report_filename}"',
                                    width=90, subsequent_indent=' ' * 11, break_on_forward_slash=True)

    logger.info(f'{fill}')

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'tree file\tNumber of tree tips\tNumber seqs in original alignment\n')

        for tree_name, stats in sorted(alignment_filtering_dict.items()):
            report_handle.write(f'{tree_name}\t{stats[0]}\t{stats[1]}\n')


def main(args,
         report_directory,
         algorithm_suffix=None,
         logger=None):
    """
    Entry point for the paragone_main.py script

    :param args: argparse namespace with subparser options for function main()
    :param str report_directory: path to directory for report files
    :param str algorithm_suffix: suffix for algorithm trees to extract fasta from
    :param logging.Logger logger: a logger object
    :return:
    """

    logger.debug(f'{"[INFO]:":10} Module fasta_from_tree was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    logger.info('')

    trimmed = False  # set default, only relevant for paralogy pruning algorithms
    if algorithm_suffix in ['mo', 'mi', 'rt']:
        if os.path.isdir('12_pre_paralog_resolution_alignments_trimmed'):
            original_alignments_directory = '12_pre_paralog_resolution_alignments_trimmed'
            trimmed = True
        else:
            original_alignments_directory = '11_pre_paralog_resolution_alignments'
            trimmed = False

    if args.from_cut_deep_paralogs:
        original_alignments_directory = args.mask_tips_alignment_directory
        treefile_directory = '08_trees_trimmed_masked_cut'
        tree_file_suffix = '.subtree'
        output_folder = f'09_sequences_from_qc_trees'

        logger.info(f'{"[INFO]:":10} ======> RECOVERING FASTA SEQUENCES CORRESPONDING TO QC TREES <======\n')

    if algorithm_suffix in ['mo']:
        treefile_directory = '14_pruned_MO'
        tree_file_suffix = '.tre'
        output_folder = f'17_selected_sequences_MO'

        logger.info(f'{"[INFO]:":10} ======> RECOVERING FASTA SEQUENCES CORRESPONDING TO MO TREES <======\n')

    if algorithm_suffix in ['mi']:
        treefile_directory = '15_pruned_MI'
        tree_file_suffix = '.tre'
        output_folder = f'18_selected_sequences_MI'

        logger.info(f'{"[INFO]:":10} ======> RECOVERING FASTA SEQUENCES CORRESPONDING TO MI TREES <======\n')

    if algorithm_suffix in ['rt']:
        treefile_directory = '16_pruned_RT'
        tree_file_suffix = '.tre'
        output_folder = f'19_selected_sequences_RT'

        logger.info(f'{"[INFO]:":10} ======> RECOVERING FASTA SEQUENCES CORRESPONDING TO RT TREES <======\n')

    # Checking input directories and files:
    directory_suffix_dict = {treefile_directory: tree_file_suffix}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    utils.createfolder(output_folder)

    # Recover alignments with sequences corresponding to tree tip names:
    filtered_alignments_folder,  alignment_filtering_dict = \
        subsample_alignments(treefile_directory,
                             output_folder,
                             tree_file_suffix,
                             original_alignments_directory,
                             trimmed=trimmed,
                             from_cut_deep_paralogs=args.from_cut_deep_paralogs,
                             logger=logger)

    # Write a report of pre-and-post filtering stats for each tree/alignments:
    write_fasta_from_tree_report(alignment_filtering_dict,
                                 report_directory,
                                 args.from_cut_deep_paralogs,
                                 algorithm_suffix,
                                 logger=logger)

    logger.info(f'{"[INFO]:":10} Finished extracting fasta sequences corresponding to tree tips.')
