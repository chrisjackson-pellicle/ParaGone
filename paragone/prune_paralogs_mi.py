#!/usr/bin/env python

# Author: # Author: Yang and Smith
# Modified by: Chris Jackson chris.jackson@rbg.vic.gov.au

"""
Input: homolog trees
Output: individual orthologs trees

If a tip is longer than the LONG_TIP_CUTOFF and also longer than 10 times its sister, cut it off. This is to fix the
leftover trees that frequently has some long tips in it.
"""

import os
import sys
import textwrap
import glob
from collections import defaultdict
import shutil

from paragone import newick3
from paragone import tree_utils
from paragone import utils
from paragone import trim_tree_tips


def write_mi_report(report_directory,
                    tree_stats_collated_dict,
                    logger=None):
    """
    Writes a *.tsv report detailing for Maximum Inclusion (MI) pruning process.

    :param str report_directory: path to directory for report files
    :param tree_stats_collated_dict: dictionary of treename:{stats}
    :param logging.Logger logger: a logger object
    :return:
    """

    report_filename = f'{report_directory}/MI_report.tsv'

    logger.info('')
    fill = utils.fill_forward_slash(f'{"[INFO]:":10} Writing Maximum Inclusion (MI) report to file: "'
                                    f'{report_filename}"',
                                    width=90, subsequent_indent=' ' * 11, break_on_forward_slash=True)

    logger.info(f'{fill}')

    trees_with_unrecognised_names_count = 0
    trees_with_fewer_than_min_ingroup_taxa_count = 0
    trees_with_1to1_orthologs_count = 0
    trees_with_clades_with_fewer_than_min_ingroup_taxa = 0
    trees_with_clades_with_greater_than_min_ingroup_taxa = 0
    trees_with_pruned_ortholog_nodes_above_trim_relative_cutoff = 0
    trees_with_pruned_ortholog_nodes_above_trim_absolute_cutoff = 0
    trees_with_pruned_ortholog_above_min_taxa = 0
    trees_with_pruned_ortholog_below_min_taxa = 0

    all_tree_stats_for_report = []

    for tree_name, dictionaries in tree_stats_collated_dict.items():

        tree_stats = [tree_name]

        try:
            check = dictionaries['unrecognised_names']
            trees_with_unrecognised_names_count += 1
            tree_stats.append(', '.join(check))
        except KeyError:
            tree_stats.append('None')

        try:
            check = dictionaries['fewer_than_min_ingroup_taxa']
            trees_with_fewer_than_min_ingroup_taxa_count += 1
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries['1to1_orthologs']
            trees_with_1to1_orthologs_count += 1
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries['clades_with_fewer_than_min_taxa']
            if len(check) != 0:
                trees_with_clades_with_fewer_than_min_ingroup_taxa += 1
                tree_stats.append(len(check))
            else:
                tree_stats.append('0')
        except KeyError:
            tree_stats.append('0')

        try:
            check = dictionaries['clades_with_greater_than_min_taxa']
            if len(check) != 0:
                trees_with_clades_with_greater_than_min_ingroup_taxa += 1
                tree_stats.append(len(check))
            else:
                tree_stats.append('0')
        except KeyError:
            tree_stats.append('0')

        try:
            check = dictionaries['pruned_ortholog_nodes_above_trim_relative_cutoff']
            tips_trimmed = False
            for ortho, tip_dict in check.items():
                if len(tip_dict) != 0:
                    tips_trimmed = True
                    break

            if tips_trimmed:
                print(check)
                trees_with_pruned_ortholog_nodes_above_trim_relative_cutoff += 1
                ortho_suffix = 1
                ortho_stats = []
                for ortho, tip_dict in check.items():
                    ortho_stats.append(f'ortho_{str(ortho_suffix)}, {len(tip_dict)} tips;')
                    ortho_suffix += 1
                tree_stats.append(' '.join(ortho_stats).rstrip(';'))

            else:
                tree_stats.append('0')
        except KeyError:
            tree_stats.append('0')

        try:
            check = dictionaries['pruned_ortholog_nodes_above_trim_absolute_cutoff']
            tips_trimmed = False
            for ortho, tip_dict in check.items():
                if len(tip_dict) != 0:
                    tips_trimmed = True
                    break

            if tips_trimmed:
                trees_with_pruned_ortholog_nodes_above_trim_absolute_cutoff += 1
                ortho_suffix = 1
                ortho_stats = []
                for ortho, tip_dict in check.items():
                    ortho_stats.append(f'ortho_{str(ortho_suffix)}, {len(tip_dict)} tips;')
                    ortho_suffix += 1
                tree_stats.append(' '.join(ortho_stats).rstrip(';'))

            else:
                tree_stats.append('0')
        except KeyError:
            tree_stats.append('0')

        try:
            check = dictionaries['pruned_ortholog_above_min_taxa']
            assert len(check) > 0
            trees_with_pruned_ortholog_above_min_taxa += 1
            tree_stats.append(len(check))
        except AssertionError:
            tree_stats.append('0')
        except KeyError:
            tree_stats.append('0')

        try:
            check = dictionaries['pruned_ortholog_below_min_taxa']
            assert len(check) > 0
            trees_with_pruned_ortholog_below_min_taxa += 1
            tree_stats.append(len(check))
        except AssertionError:
            tree_stats.append('0')
        except KeyError:
            tree_stats.append('0')

        all_tree_stats_for_report.append(tree_stats)

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'\t'
                            f'Unrecognised taxa (tree skipped)\t'
                            f'< than minimum ingroup taxa (tree skipped)\t'
                            f'1-to-1 orthologs\t'
                            f'Ortholog clades < minimum ingroup taxa\t'
                            f'Ortholog clades > minimum ingroup taxa\t'
                            f'Pruned orthologs with tips > relative cutoff\t'
                            f'Pruned orthologs with tips > absolute cutoff\t'
                            f'Pruned orthologs > than minimum taxa\t'
                            f'Pruned orthologs < than minimum taxa'
                            f'\n')

        report_handle.write(f'Number of trees\t'
                            f'{trees_with_unrecognised_names_count}\t'
                            f'{trees_with_fewer_than_min_ingroup_taxa_count}\t'
                            f'{trees_with_1to1_orthologs_count}\t'
                            f'{trees_with_clades_with_fewer_than_min_ingroup_taxa}\t'
                            f'{trees_with_clades_with_greater_than_min_ingroup_taxa}\t'
                            f'{trees_with_pruned_ortholog_nodes_above_trim_relative_cutoff}\t'
                            f'{trees_with_pruned_ortholog_nodes_above_trim_absolute_cutoff}\t'
                            f'{trees_with_pruned_ortholog_above_min_taxa}\t'
                            f'{trees_with_pruned_ortholog_below_min_taxa}'
                            f'\n')

        for stats in all_tree_stats_for_report:
            stats_joined = '\t'.join([str(stat) for stat in stats])
            report_handle.write(f'{stats_joined}\n')


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

    logger.debug(f'{"[INFO]:":10} Module prune_paralogs_mi was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]),
                         width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11, break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    logger.info('')
    logger.info(f'{"[INFO]:":10} ======> PRUNING PARALOGS WITH MI ALGORITHM <======\n')

    # Checking input directories and files:
    treefile_directory = '13_pre_paralog_resolution_trees'
    tree_file_suffix = '.treefile'
    directory_suffix_dict = {treefile_directory: tree_file_suffix}
    in_and_outgroups_list = '00_logs_and_reports/reports/in_and_outgroups_list.tsv'
    file_list = [in_and_outgroups_list]

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Create output folder for pruned trees:
    output_folder = f'15_pruned_MI'
    utils.createfolder(output_folder)

    # Parse the ingroup and outgroup text file:
    ingroups, outgroups = utils.parse_ingroup_and_outgroup_file(in_and_outgroups_list,
                                                                logger=logger)

    # Create dict for report file:
    tree_stats_collated = defaultdict(lambda: defaultdict())

    # Iterate over tree and prune with MO algorithm:
    for treefile in glob.glob(f'{treefile_directory}/*{tree_file_suffix}'):
        treefile_basename = os.path.basename(treefile)
        output_file_id = f'{output_folder}/{tree_utils.get_cluster_id(treefile_basename)}'

        logger.info(f'{"[INFO]:":10} Analysing tree {treefile_basename}...')

        with open(treefile, "r") as infile:
            intree = newick3.parse(infile.readline())
            curroot = intree
            names = tree_utils.get_front_names(curroot)
            num_tips, num_taxa = len(names), len(set(names))
            ingroup_names = []
            outgroup_names = []
            unrecognised_names = []

            for name in names:
                if name in ingroups:
                    ingroup_names.append(name)
                elif name in outgroups:
                    outgroup_names.append(name)
                else:
                    unrecognised_names.append(name)

            # Check for unrecognised tip names and skip tree if present:
            if unrecognised_names:
                fill = textwrap.fill(
                    f'{"[WARNING]:":10} Taxon names {unrecognised_names} in tree {treefile_basename} not found in '
                    f'ingroups or outgroups. Skipping tree...',
                    width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

                logger.warning(f'{fill}')

                tree_stats_collated[treefile_basename]['unrecognised_names'] = unrecognised_names
                continue

            # Check if tree contains more than the minimum number of taxa:
            if len(ingroup_names) < args.minimum_taxa:
                fill = textwrap.fill(
                    f'{"[WARNING]:":10} Tree {treefile_basename} contains {len(ingroup_names)} ingroup taxa; '
                    f'minimum_taxa required is {args.minimum_taxa}. Skipping tree...',
                    width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

                logger.warning(f'{fill}')

                tree_stats_collated[treefile_basename]['fewer_than_min_ingroup_taxa'] = newick3.tostring(curroot)
                continue

            # Check if paralogs are present; if not, write 1to1ortho tree (optional):
            if tree_utils.get_front_score(curroot) >= args.minimum_taxa:  # get_front_score -1 if paralogs present
                fill = textwrap.fill(
                    f'{"[INFO]:":10} Tree {treefile_basename} contain no duplicated taxon names (i.e. paralogs).',
                    width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

                logger.info(f'{fill}')

                tree_stats_collated[treefile_basename]['1to1_orthologs'] = newick3.tostring(curroot)

                if not args.ignore_1to1_orthologs:

                    fill = textwrap.fill(
                        f'{"[INFO]:":10} Writing tree {treefile_basename} to {output_file_id}.1to1ortho.tre',
                        width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

                    logger.info(f'{fill}')

                    shutil.copy(treefile, f'{output_file_id}.1to1ortho.tre')
                else:
                    fill = textwrap.fill(
                        f'{"[INFO]:":10} Parameter --ignore_1to1_orthologs provided. Skipping tree '
                        f'{treefile_basename}...',
                        width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

                    logger.info(f'{fill}')

                continue

            # If duplicate taxon names present, prune tree:
            going = True
            pruned_clades = []
            clades_with_fewer_than_min_taxa = []
            clades_with_greater_than_min_taxa = []

            while going:
                highest = 0
                highest_node = None
                scores_dict = {}  # key is node, value is a tuple (front_score,back_score)

                for node in curroot.iternodes():
                    # front_score and back_score (below) are > 0 if no paralogs in current clade, or -1 if paralogs
                    front_score = int(tree_utils.get_front_score(node))
                    back_score = int(tree_utils.get_back_score(node, curroot))
                    scores_dict[node] = (front_score, back_score)

                    if front_score > highest or back_score > highest:
                        highest_node = node  # node with the greatest number of non-duplicated leaf names
                        highest = max(front_score, back_score)  # number of taxa in the clade with the greatest
                        # number of non-duplicated leaf names.

                # We've now identified the node in the tree with the greatest number of non-repeating taxa in either
                # the front or back clade, but this number includes any outgroup sequences. We need to determine which
                # clade (front or back) and count ingroup taxa only:
                highest_front, highest_back = scores_dict[highest_node]

                if highest_front > highest_back:
                    highest_ingroup_names_mi = tree_utils.get_front_ingroup_names(highest_node, ingroups)
                else:
                    highest_ingroup_names_mi = tree_utils.get_back_ingroup_names(highest_node, curroot, ingroups)

                if highest >= args.minimum_taxa:  # prune

                    if len(highest_ingroup_names_mi) >= args.minimum_taxa:
                        clades_with_greater_than_min_taxa.append(highest_node)
                    else:
                        clades_with_fewer_than_min_taxa.append(highest_node)

                    curroot, done = tree_utils.prune(scores_dict[highest_node],
                                                     highest_node,
                                                     curroot,
                                                     pruned_clades,
                                                     highest_ingroup_names_mi,
                                                     args.minimum_taxa,
                                                     logger=logger)
                    if done:
                        going = False
                        break

                    elif len(tree_utils.get_front_ingroup_names(curroot, ingroups)) < args.minimum_taxa:
                        clades_with_fewer_than_min_taxa.append(newick3.tostring(curroot))
                        going = False
                        break
                else:
                    clades_with_fewer_than_min_taxa.append(newick3.tostring(curroot))
                    going = False
                    break

            tree_stats_collated[treefile_basename]['clades_with_fewer_than_min_taxa'] = \
                clades_with_fewer_than_min_taxa

            tree_stats_collated[treefile_basename]['clades_with_greater_than_min_taxa'] = \
                clades_with_greater_than_min_taxa

            if len(pruned_clades) > 0:

                fill = utils.fill_forward_slash(f'{"[INFO]:":10} {len(pruned_clades)} pruned clades recovered for '
                                                f'tree {treefile_basename}',
                                                width=90, subsequent_indent=' ' * 11, break_on_forward_slash=True)
                logger.info(f'{fill}')

                count = 1

                pruned_orthologs_relative_cutoffs = {}
                pruned_orthologs_absolute_cutoffs = {}
                pruned_orthologs_above_minimum_taxa = []
                pruned_orthologs_below_minimum_taxa = []

                for tree in pruned_clades:
                    if tree.nchildren == 2:
                        node, tree = tree_utils.remove_kink(tree, tree)

                    # Trim tips:
                    tree, nodes_above_absolute_cutoff, nodes_above_relative_cutoff = \
                        trim_tree_tips.trim(tree,
                                            args.relative_tip_cutoff,
                                            args.absolute_tip_cutoff,
                                            args.minimum_taxa,
                                            tree_name=treefile_basename,
                                            logger=logger)

                    pruned_orthologs_relative_cutoffs[tree] = nodes_above_relative_cutoff
                    pruned_orthologs_absolute_cutoffs[tree] = nodes_above_absolute_cutoff

                    # Write pruned ortholog trees if above minimum taxa:
                    if tree:
                        ingroup_names_mi = tree_utils.get_front_ingroup_names(tree, ingroups)
                        logger.debug(f'Ingroup taxa in ortho after MI pruning: {ingroup_names_mi}')

                        if len(ingroup_names_mi) >= args.minimum_taxa:
                            tree_output_filename = f'{output_file_id}.MIortho{str(count)}.tre'

                            pruned_orthologs_above_minimum_taxa.append(tree)

                            with open(tree_output_filename, "w") as outfile:
                                outfile.write(newick3.tostring(tree) + ";\n")

                            count += 1
                        else:
                            logger.info(f'{"[INFO]:":10} After trimming tips, pruned ortholog from {treefile_basename} '
                                        f'contained fewer than minimum ingroup taxa value of {args.minimum_taxa}. '
                                        f'Skipping pruned ortholog...')

                            pruned_orthologs_below_minimum_taxa.append(tree)

                # Recover stats in dictionary for report:
                tree_stats_collated[treefile_basename]['pruned_ortholog_nodes_above_trim_relative_cutoff'] = \
                    pruned_orthologs_relative_cutoffs

                tree_stats_collated[treefile_basename]['pruned_ortholog_nodes_above_trim_absolute_cutoff'] = \
                    pruned_orthologs_absolute_cutoffs

                tree_stats_collated[treefile_basename]['pruned_ortholog_above_min_taxa'] = \
                    pruned_orthologs_above_minimum_taxa

                tree_stats_collated[treefile_basename]['pruned_ortholog_below_min_taxa'] = \
                    pruned_orthologs_below_minimum_taxa

    # Write a *.tsv report file:
    write_mi_report(report_directory,
                    tree_stats_collated,
                    logger=logger)

    fill = textwrap.fill(f'{"[INFO]:":10} Finished extracting putative ortholog trees using the Maximum Inclusion (MI) '
                         f'algorithm. Trees have been written to directory: "{output_folder}".',
                         width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

    logger.info(f'{fill}')