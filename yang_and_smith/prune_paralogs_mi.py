#!/usr/bin/env python

# Author: # Author: Yang and Smith
# Modified by: Chris Jackson chris.jackson@rbg.vic.gov.au

"""
Input: homolog trees
Output: individual orthologs trees

if a tip is longer than the LONG_TIP_CUTOFF
and also long than 10 times its sister, cut it off
This is to fix the leftover trees that frequently has some long tips in it

If not to output 1-to-1 orthologs, for example, already analysed these
set OUTPUT_1to1_ORTHOLOGS to False
"""

import os
import sys
import textwrap
import glob
from collections import defaultdict
import shutil

from yang_and_smith import newick3
from yang_and_smith import tree_utils
from yang_and_smith import utils
from yang_and_smith import trim_tree_tips


def write_mi_report(treefile_directory,
                    tree_stats_collated_dict,
                    logger=None):
    """
    Writes a *.tsv report detailing for Maximum Inclusion (MI) pruning process.

    :param str treefile_directory: name of tree file directory for report filename
    :param tree_stats_collated_dict: dictionary of treename:{stats}
    :param logging.Logger logger: a logger object
    :return:
    """

    basename = os.path.basename(treefile_directory)
    report_filename = f'{basename}_MI_report.tsv'

    logger.info(f'{"[INFO]:":10} Writing Maximum Inclusion (MI) report to file {report_filename}')

    trees_with_unrecognised_names_count = 0
    trees_with_fewer_than_min_ingroup_taxa_count = 0
    trees_with_1to1_orthologs_count = 0
    trees_with_no_outgroup_taxa_count = 0
    tree_with_duplicate_taxa_in_outgroup_count = 0
    trees_with_monophyletic_outgroups_count = 0
    trees_with_non_monophyletic_outgroups_count = 0
    trees_with_mo_output_file_above_minimum_taxa_count = 0
    trees_with_mo_output_file_below_minimum_taxa_count = 0

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
            check = dictionaries['no_outgroup_taxa']
            trees_with_no_outgroup_taxa_count += 1
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries['duplicate_taxa_in_outgroup']
            tree_with_duplicate_taxa_in_outgroup_count += 1
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries['monophyletic_outgroups']
            trees_with_monophyletic_outgroups_count += 1
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries['non_monophyletic_outgroups']
            trees_with_non_monophyletic_outgroups_count += 1
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries['mo_output_file_above_minimum_taxa']
            trees_with_mo_output_file_above_minimum_taxa_count += 1
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries['mo_output_file_below_minimum_taxa']
            trees_with_mo_output_file_below_minimum_taxa_count += 1
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        all_tree_stats_for_report.append(tree_stats)

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'\t'
                            f'Unrecognised taxa (tree skipped)\t'
                            f'< than minimum ingroup taxa (tree skipped)\t'
                            f'1-to-1 orthologs\t'
                            f'No outgroup taxa\t'
                            f'Duplicate taxa in the outgroup\t'
                            f'Pputative paralogs and monophyletic outgroup\t'
                            f'Putative paralogs and non-monophyletic outgroup\t'
                            f'MO pruned trees > than minimum taxa\t'
                            f'MO pruned trees < than minimum taxa'
                            f'\n')

        report_handle.write(f'Number of trees\t'
                            f'{trees_with_unrecognised_names_count}\t'
                            f'{trees_with_fewer_than_min_ingroup_taxa_count}\t'
                            f'{trees_with_1to1_orthologs_count}\t'
                            f'{trees_with_no_outgroup_taxa_count}\t'
                            f'{tree_with_duplicate_taxa_in_outgroup_count}\t'
                            f'{trees_with_monophyletic_outgroups_count}\t'
                            f'{trees_with_non_monophyletic_outgroups_count}\t'
                            f'{trees_with_mo_output_file_above_minimum_taxa_count}\t'
                            f'{trees_with_mo_output_file_below_minimum_taxa_count}'
                            f'\n')

        for stats in all_tree_stats_for_report:
            stats_joined = '\t'.join([str(stat) for stat in stats])
            report_handle.write(f'{stats_joined}\n')


def main(args):
    """
    Entry point for the resolve_paralogs.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = utils.setup_logger(__name__, 'logs_resolve_paralogs/11_prune_paralogs_MI')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand prune_paralogs_mi was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')
    logger.debug(args)

    # Create output folder for pruned trees:
    output_folder = f'{os.path.basename(args.treefile_directory)}_pruned_MI'
    utils.createfolder(output_folder)

    # Create dict for report file:
    tree_stats_collated = defaultdict(lambda: defaultdict())

    # Iterate over tree and prune with MO algorithm:
    for treefile in glob.glob(f'{args.treefile_directory}/*{args.tree_file_suffix}'):
        treefile_basename = os.path.basename(treefile)
        output_file_id = f'{output_folder}/{tree_utils.get_cluster_id(treefile_basename)}'

        logger.info(f'{"[INFO]:":10} Analysing tree {treefile_basename}...')

        with open(treefile, "r") as infile:
            intree = newick3.parse(infile.readline())
            curroot = intree

            pruned_clades = []

            # Check if paralogs are present; if not, write 1to1ortho tree (optional):
            if tree_utils.get_front_score(curroot) >= args.minimum_taxa:
                logger.info(f'{"[INFO]:":10} Tree {treefile_basename} contain no duplicated taxon names (i.e. '
                            f'paralogs).')

                tree_stats_collated[treefile_basename]['1to1_orthologs'] = newick3.tostring(curroot)

                if not args.ignore_1to1_orthologs:
                    logger.info(f'{"[INFO]:":10} Writing tree {treefile_basename} to {output_file_id}.1to1ortho.tre')
                    shutil.copy(treefile, f'{output_file_id}.1to1ortho.tre')
                else:
                    logger.info(f'{"[INFO]:":10} Parameter --ignore_1to1_orthologs provided. Skipping tree...'
                                f' {treefile_basename}')

            # If duplicate taxon names present, prune tree:
            else:
                going = True
                pruned_clades = []
                clades_with_fewer_than_min_taxa = []

                while going:
                    highest = 0
                    highest_node = None
                    scores_dict = {}  # key is node, value is a tuple (front_score,back_score)

                    for node in curroot.iternodes():
                        # front_score and back_score below > 0 if no paralogs in current clade, or -1 if paralogs
                        front_score = int(tree_utils.get_front_score(node))
                        back_score = int(tree_utils.get_back_score(node, curroot))
                        scores_dict[node] = (front_score, back_score)

                        if front_score > highest or back_score > highest:
                            highest_node = node  # node with the greatest number of non-duplicated leaf names
                            highest = max(front_score, back_score)  # number of taxa in the clade with the greatest
                            # number of non-duplicated leaf names.

                    if highest >= args.minimum_taxa:  # prune

                        curroot, done = tree_utils.prune(scores_dict[highest_node],
                                                         highest_node,
                                                         curroot,
                                                         pruned_clades,
                                                         logger=logger)

                        if done:
                            going = False
                            break
                        elif len(curroot.leaves()) < args.minimum_taxa:
                            print(newick3.tostring(curroot))
                            clades_with_fewer_than_min_taxa.append(newick3.tostring(curroot))
                            going = False
                            break
                    else:
                        going = False
                        break

                tree_stats_collated[treefile_basename]['clades_with_fewer_than_min_taxa'] = \
                    clades_with_fewer_than_min_taxa

        if len(pruned_clades) > 0:
            logger.info(f'{"[INFO]:":10} {len(pruned_clades)} pruned clades recovered for tree {treefile_basename}')

            tree_stats_collated[treefile_basename]['trees_with_pruned_clades'] = pruned_clades
            count = 1

            pruned_orthologs_absolute_and_relative_cutoffs = {}
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
                                        tree_name=treefile_basename,
                                        logger=logger)

                pruned_orthologs_absolute_and_relative_cutoffs[tree] = \
                    (nodes_above_absolute_cutoff, nodes_above_relative_cutoff)

                # Write pruned ortholog trees if above minimum taxa:
                if tree:
                    if len(tree.leaves()) >= args.minimum_taxa:
                        tree_output_filename = f'{output_file_id}.MIortho{str(count)}.tre'

                        pruned_orthologs_above_minimum_taxa.append(tree)

                        with open(tree_output_filename, "w") as outfile:
                            outfile.write(newick3.tostring(tree) + ";\n")

                        count += 1
                    else:
                        logger.info(f'{"[INFO]:":10} After trimming tips, pruned ortholog from {treefile_basename} '
                                    f'contained fewer than minimum taxa value of {args.minimum_taxa}. Skipping pruned '
                                    f'ortholog...')

                        pruned_orthologs_below_minimum_taxa.append(tree)

            tree_stats_collated[treefile_basename]['pruned_orthologs_nodes_above_trim_cutoffs'] = \
                pruned_orthologs_absolute_and_relative_cutoffs

            tree_stats_collated[treefile_basename]['pruned_orthologs_nodes_above_min_taxa'] = \
                pruned_orthologs_above_minimum_taxa

            tree_stats_collated[treefile_basename]['pruned_orthologs_nodes_below_min_taxa'] = \
                pruned_orthologs_below_minimum_taxa

    # Write a *.tsv report file:
    write_mi_report(args.treefile_directory,
                    tree_stats_collated,
                    logger=logger)

