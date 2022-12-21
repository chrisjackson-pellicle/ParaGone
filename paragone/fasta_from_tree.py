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
                         from_cut_deep_paralogs=False,
                         algorithm_suffix=None,
                         logger=None):
    """
    Takes a pruned/QC'd tree file, finds the original matching alignment, and sub-samples that alignment to recover
    only sequences corresponding to tree tip names.

    :param str treefile_directory: path to directory containing tree newick files
    :param str output_folder: path to output folder for selected fasta alignments
    :param str tree_suffix: suffix for the tree files
    :param str alignment_directory: path to the directory containing fasta alignments
    :param bool from_cut_deep_paralogs: if True, process tree file names accordingly to recover gene names
    :param None/str algorithm_suffix: if extracting seqs from pruned trees, the algorithm suffix mo/rt/mi, else None
    :param logging.Logger logger: a logger object
    :return str, dict output_folder, alignment_filtering_dict: path the output folder with filtered alignments,
    dictionary of filtering stats for each tree/alignment
    """

    logger.info(f'{"[INFO]:":10} Recovering alignment sequences corresponding to tree tip names...')

    # Check for trim/clean status of original alignments:
    if re.search('hmmcleaned', alignment_directory):
        logger.debug(f'Input alignment folder is {alignment_directory}; sequenced were trimmed and cleaned')
        alignment_suffix = '.aln.trimmed.hmm.fasta'
    elif re.search('trimmed', alignment_directory):
        logger.debug(f'Input alignment folder is {alignment_directory}; sequenced were trimmed but not cleaned')
        alignment_suffix = '.aln.trimmed.fasta'
    else:
        logger.debug(f'Input alignment folder is {alignment_directory}; sequenced were neither trimmed or cleaned')
        alignment_suffix = '.aln.fasta'

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
            matching_alignment = f'{alignment_directory}/{alignment_prefix}.outgroup_added.aln.trimmed.fasta'

        # Read in original alignments and select seqs matching tree termini:
        alignment = AlignIO.read(matching_alignment, "fasta")
        subalignment = Bio.Align.MultipleSeqAlignment([])
        for k in range(0, len(alignment)):
            for j in range(0, len(tree_terminals)):
                if tree_terminals[j].name == alignment[k].id:
                    subalignment.append(alignment[k])

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
    Writes a *.tsv report detailing number of tips in QC'd tree, number of sequences in original QC'd alignment,
    and number of sequences in filtered alignment.

    :param dict alignment_filtering_dict: dictionary of filtering stats for each tree/alignment
    :param str report_directory: path to directory for report files
    :param bool from_cut_deep_paralogs: if True, add 'cut' to report filename
    :param strNone algorithm_suffix: if extracting seqs from pruned trees, the algorithm suffix mo/rt/mi, else None
    :param logging.Logger logger: a logger object
    :return:
    """

    if from_cut_deep_paralogs:
        report_filename = f'{report_directory}/fasta_from_qc_trees_report.tsv'
    else:
        # report_filename = f'00_logs_and_reports_resolve_paralogs/reports/' \
        #                   f'{re.sub("^[0-9]{2}_", "", basename)}_fasta_from_tree_{algorithm_suffix}_report.tsv'
        report_filename = f'00_logs_and_reports_resolve_paralogs/reports/fasta_from_tree_{algorithm_suffix}_report.tsv'

    fill = utils.fill_forward_slash(f'{"[INFO]:":10} Writing trim tips report to file: "{report_filename}"',
                                    width=90, subsequent_indent=' ' * 11, break_on_forward_slash=True)

    logger.info(f'{fill}')

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'QC tree file\tNumber of tree tips\tNumber seqs in original alignment\tNumber seqs '
                            f'filtered alignment\n')

        for tree_name, stats in alignment_filtering_dict.items():
            report_handle.write(f'{tree_name}\t{stats[0]}\t{stats[1]}\t{stats[2]}\n')


def main(args,
         report_directory,
         logger=None):
    """
    Entry point for the paragone_main.py script

    :param args: argparse namespace with subparser options for function main()
    :param str report_directory: path to directory for report files
    :param logging.Logger logger: a logger object
    :return:
    """

    logger.debug(f'{"[INFO]:":10} Module fasta_from_tree was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    logger.info('')
    logger.info(f'{"[INFO]:":10} ======> RECOVERING FASTA SEQUENCES CORRESPONDING TO QC TREES <======\n')

    algorithm_suffix = None

    # Checking input directories and files:
    directory_suffix_dict = {args.treefile_directory: args.tree_file_suffix}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Create output folder:
    if args.from_cut_deep_paralogs:
        output_folder = f'09_alignments_from_qc_trees'
    elif args.from_prune_paralogs:
        output_folder = f'22_selected_sequences_{algorithm_suffix}_batches'

    utils.createfolder(output_folder)

    # Recover alignments with sequences corresponding to tree tip names:
    filtered_alignments_folder,  alignment_filtering_dict = \
        subsample_alignments(args.treefile_directory,
                             output_folder,
                             args.tree_file_suffix,
                             args.mask_tips_alignment_directory,  # same alignments used for mask_tips step
                             from_cut_deep_paralogs=args.from_cut_deep_paralogs,
                             algorithm_suffix=algorithm_suffix,
                             logger=logger)

    # Write a report of pre-and-post filtering stats for each tree/alignments:
    write_fasta_from_tree_report(alignment_filtering_dict,
                                 report_directory,
                                 args.from_cut_deep_paralogs,
                                 algorithm_suffix,
                                 logger=logger)

    logger.info(f'{"[INFO]:":10} Finished extracting fasta sequences corresponding to tree tips.')
