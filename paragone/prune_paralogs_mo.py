#!/usr/bin/env python

# Author: Yang and Smith (2014), modified by Alexander Schmidt-Lebuhn

# Modified by: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
Taxon duplication? --No--> output one-to-one orthologs
        |
       Yes
        |
Outgroup present? --No--> ignore this homolog
        |
       Yes
        |
Outgroup taxon duplication? --Yes--> ignore this homolog
        |
        No
        |
Outgroup monophyletic? --No--> ignore this homolog
        |
       Yes
        |
Infer orthologs by using monophyletic, non-repeating outgroups

If not to output 1-to-1 orthologs, for example, already analysed these
set OUTPUT_1to1_ORTHOLOGS to False
"""

import os
import sys
from collections import defaultdict
import glob
import shutil
import textwrap

from paragone import utils
from paragone import tree_utils
from paragone import phylo3
from paragone import newick3


def reroot_with_monophyletic_outgroups(root,
                                       outgroups,
                                       logger=None):
    """
    Check if outgroups are monophyletic and non-repeating and reroot, otherwise return None

    :param phylo3.Node root: tree object parsed by newick3.parse
    :param list outgroups: a list of outgroup taxon names present in the tree
    :param logging.Logger logger: a logger object
    :return:
    """

    leaves = root.leaves()
    outgroup_matches = {}  # Key is label, value is the tip node object

    # Since no taxon repeat in outgroups name and leaf is one-to-one  # CJJ what the latter clause mean?
    outgroup_labels = []
    for leaf in leaves:
        label = leaf.label  # e.g. 376678.main or 376728.0, etc
        name = tree_utils.get_name(label)  # e.g. 376678 or 376728, etc
        if name in outgroups:
            outgroup_matches[label] = leaf
            outgroup_labels.append(label)

    if len(outgroup_labels) == 1:  # Tree contains a single outgroup sequence
        # cannot reroot on a tip so have to go one more node into the ingroup:
        new_root = outgroup_matches[outgroup_labels[0]].parent
        return phylo3.reroot(root, new_root)

    else:  # Tree has multiple outgroup sequences. Check monophyly and reroot:
        newroot = None
        for node in root.iternodes():  # Iterate over nodes and try to find one with monophyletic outgroup
            if node == root:
                continue  # Skip the root

            front_names = tree_utils.get_front_names(node)
            back_names = tree_utils.get_back_names(node, root)
            front_in_names, front_out_names, back_in_names, back_out_names = 0, 0, 0, 0

            # Get counts of ingroup and outgroup taxa at front and back of the current node:
            for i in front_names:
                if i in outgroups:
                    front_out_names += 1
                else:
                    front_in_names += 1
            for j in back_names:
                if j in outgroups:
                    back_out_names += 1
                else:
                    back_in_names += 1

            if front_in_names == 0 and front_out_names > 0 and back_in_names > 0 and back_out_names == 0:
                newroot = node.parent  # ingroup at back, outgroup in front CJJ added.parent - bugfix?
                break

            if front_in_names > 0 and front_out_names == 0 and back_in_names == 0 and back_out_names > 0:
                newroot = node.parent  # ingroup in front, outgroup at back
                break

        if newroot:
            return phylo3.reroot(root, newroot)
        else:
            return None


def prune_paralogs_from_rerooted_homotree(root,
                                          outgroups,
                                          logger=None):
    """
    Prunes a tree containing monophletic outgroup sequences to recover the ingroup clade with the largest number of
    non-repeating taxon names. Returns a tree containing the outgroup sequences as well as ingroup sequences.

    :param phylo3.Node root: tree object parsed by newick3.parse
    :param list outgroups: list of outgroup names recovered from in_and_outgroup_list file
    :param logging.Logger logger: a logger object
    :return phylo3.Node root: tree object after pruning with Monophyletic Outgroups (MO) algorithm
    """

    if len(tree_utils.get_front_names(root)) == len(set(tree_utils.get_front_names(root))):
        return root  # no pruning needed CJJ This is same as 1to1_orthologs, isn't it?

    # Check for duplications at the root first. One or two of the trifurcating root clades are ingroup clades:
    node0, node1, node2 = root.children[0], root.children[1], root.children[2]
    out0, out1, out2 = len(tree_utils.get_front_outgroup_names(node0, outgroups)),\
                       len(tree_utils.get_front_outgroup_names(node1, outgroups)),\
                       len(tree_utils.get_front_outgroup_names(node2, outgroups))

    logger.debug(f'Outgroup taxon count in node0, node1, node2 is: {out0}, {out1}, {out2}')

    # Identify the ingroup clades and check for names overlap:
    if out0 == 0 and out1 == 0:  # 0 and 1 are the ingroup clades
        name_set0 = set(tree_utils.get_front_names(node0))
        name_set1 = set(tree_utils.get_front_names(node1))
        if len(name_set0.intersection(name_set1)) > 0:

            if len(name_set0) > len(name_set1):  # cut the side with fewer taxa
                logger.debug(f'Cutting node1: {newick3.tostring(node1)}')
                root.remove_child(node1)
                node1.prune()
            else:
                root.remove_child(node0)  # CJJ arbitrary removal of node0 rather than node1 if same number taxa?
                logger.debug(f'Cutting node0: {newick3.tostring(node0)}')
                node0.prune()

    elif out1 == 0 and out2 == 0:  # 1 and 2 are the ingroup clades
        name_set1 = set(tree_utils.get_front_names(node1))
        name_set2 = set(tree_utils.get_front_names(node2))
        if len(name_set1.intersection(name_set2)) > 0:
            if len(name_set1) > len(name_set2):  # cut the side with fewer taxa
                logger.debug(f'Cutting node2: {newick3.tostring(node2)}')
                root.remove_child(node2)
                node2.prune()
            else:
                root.remove_child(node1)
                logger.debug(f'Cutting node1: {newick3.tostring(node1)}')
                node1.prune()

    elif out0 == 0 and out2 == 0:  # 0 and 2 are the ingroup clades
        name_set0 = set(tree_utils.get_front_names(node0))
        name_set2 = set(tree_utils.get_front_names(node2))
        if len(name_set0.intersection(name_set2)) > 0:
            if len(name_set0) > len(name_set2):  # cut the side with fewer taxa
                root.remove_child(node2)
                logger.debug(f'Cutting node2: {newick3.tostring(node2)}')
                node2.prune()
            else:
                root.remove_child(node0)
                logger.debug(f'Cutting node0: {newick3.tostring(node0)}')
                node0.prune()

    else:
        raise ValueError('More than one clade with outgroup sequences!')

    # If there are still taxon duplications (putative paralogs) in the ingroup clade, keep pruning:
    while len(tree_utils.get_front_names(root)) > len(set(tree_utils.get_front_names(root))):
        for node in root.iternodes(order=0):  # PREORDER, root to tip  CJJ: this tree includes outgroup taxa

            if node.istip:
                continue
            elif node == root:
                continue

            child0, child1 = node.children[0], node.children[1]
            name_set0 = set(tree_utils.get_front_names(child0))
            name_set1 = set(tree_utils.get_front_names(child1))
            if len(name_set0.intersection(name_set1)) > 0:
                if len(name_set0) > len(name_set1):  # cut the side with fewer taxa
                    node.remove_child(child1)
                    child1.prune()
                else:
                    node.remove_child(child0)
                    child0.prune()
                node, root = tree_utils.remove_kink(node, root)  # no re-rooting here
                break

    return root


def write_mo_report(report_directory,
                    tree_stats_collated_dict,
                    logger=None):
    """
    Writes a *.tsv report detailing for Monophyletic Outgroup (MO) pruning process.

    :param str treefile_directory: name of tree file directory for report filename
    :param tree_stats_collated_dict: dictionary of treename:{stats}
    :param logging.Logger logger: a logger object
    :return:
    """

    report_filename = f'{report_directory}/MO_report.tsv'

    logger.info('')
    fill = utils.fill_forward_slash(f'{"[INFO]:":10} Writing trim tips report to file: "{report_filename}"',
                                    width=90, subsequent_indent=' ' * 11, break_on_forward_slash=True)

    logger.info(f'{fill}')

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
                            f'Putative paralogs and monophyletic outgroup\t'
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

    logger.debug(f'{"[INFO]:":10} Module prune_paralogs_mo was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]),
                         width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11, break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    logger.info('')
    logger.info(f'{"[INFO]:":10} ======> PRUNING PARALOGS WITH MO ALGORITHM <======\n')

    # Checking input directories and files:
    treefile_directory = '13_pre_paralog_resolution_trees'
    tree_file_suffix = '.treefile'
    directory_suffix_dict = {treefile_directory: tree_file_suffix}
    in_and_outgroups_list = 'in_and_outgroups_list.txt'
    file_list = [in_and_outgroups_list]

    # if args.in_and_outgroup_list:
    #     file_list.append('in_and_outgroups_list.txt')

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Create output folder for pruned trees:
    output_folder = f'14_pruned_MO'
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

        # Read in the tree and check number of taxa:
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

        # If the tree has no taxon duplication, no cutting is needed:
        if num_tips == num_taxa:
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
                    f'{"[INFO]:":10} Parameter --ignore_1to1_orthologs provided. Skipping tree {treefile_basename}...',
                    width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

                logger.info(f'{fill}')

        else:
            # Now need to deal with taxon duplications. Check to make sure that the ingroup and outgroup names were
            # set correctly:
            logger.info(f'{"[INFO]:":10} Tree {treefile_basename} contains paralogs...')

            # If no outgroup at all, do not attempt to resolve paralogs:
            if len(outgroup_names) == 0:
                logger.info(f'{"[WARNING]:":10} Tree {treefile_basename} contains no outgroup taxa. Skipping tree...')

                tree_stats_collated[treefile_basename]['no_outgroup_taxa'] = newick3.tostring(curroot)

            # Skip the tree if there are duplicated outgroup taxa
            elif len(outgroup_names) > len(set(outgroup_names)):

                fill = textwrap.fill(
                    f'{"[WARNING]:":10} Tree {treefile_basename} contains duplicate taxon names in the outgroup taxa. '
                    f'Skipping tree...',
                    width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

                logger.warning(f'{fill}')

                tree_stats_collated[treefile_basename]['duplicate_taxa_in_outgroup'] = newick3.tostring(curroot)

            else:  # At least one outgroup present and no outgroup duplication
                if curroot.nchildren == 2:  # need to reroot
                    temp, curroot = tree_utils.remove_kink(curroot, curroot)

                # Check if the outgroup sequences are monophyletic:
                curroot = reroot_with_monophyletic_outgroups(curroot,
                                                             outgroups,
                                                             logger=logger)

                # Only return one tree after pruning:
                if curroot:  # i.e. the outgroup was monophyletic
                    fill = textwrap.fill(
                        f'{"[INFO]:":10} Outgroup sequences are monophyletic for tree {treefile_basename}.',
                        width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

                    logger.info(f'{fill}')

                    tree_stats_collated[treefile_basename]['monophyletic_outgroups'] = newick3.tostring(curroot)

                    # Write re-rooted trees with monophyletic outgroup to file:
                    with open(f'{output_file_id}.reroot', "w") as outfile:
                        outfile.write(newick3.tostring(curroot) + ";\n")

                    # Prune the re-rooted tree with the MO algorith:
                    fill = textwrap.fill(
                        f'{"[INFO]:":10} Applying Monophyletic Outgroup algorithm to tree {treefile_basename}...',
                        width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

                    logger.info(f'{fill}')

                    ortho = prune_paralogs_from_rerooted_homotree(curroot,
                                                                  outgroups,
                                                                  logger=logger)

                    # Filter out pruned trees that have fewer ingroup taxa than the minimum_taxa value:
                    ingroup_names_mo = tree_utils.get_front_ingroup_names(curroot, ingroups)
                    logger.debug(f'Ingroup taxa in ortho after MO pruning: {ingroup_names_mo}')

                    if len(set(ingroup_names_mo)) >= args.minimum_taxa:
                        with open(f'{output_file_id}.ortho.tre', "w") as outfile:
                            outfile.write(newick3.tostring(ortho) + ";\n")

                            tree_stats_collated[treefile_basename]['mo_output_file_above_minimum_taxa'] = \
                                newick3.tostring(curroot)
                    else:
                        fill = textwrap.fill(
                            f'{"[WARNING]:":10} After pruning with MO algorith, tree {treefile_basename} contains'
                            f' {len(set(tree_utils.get_front_names(curroot)))} taxa; parameter --minimum_taxa is'
                            f' {args.minimum_taxa}. No tree file will be written.',
                            width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

                        logger.warning(f'{fill}')

                        tree_stats_collated[treefile_basename]['mo_output_file_below_minimum_taxa'] = \
                            newick3.tostring(curroot)
                else:
                    logger.info(f'{"[INFO]:":10} Outgroup non-monophyletic for tree {treefile_basename}')

                    tree_stats_collated[treefile_basename]['non_monophyletic_outgroups'] = \
                        newick3.tostring(curroot)

    # Write a *.tsv report file:
    write_mo_report(report_directory,
                    tree_stats_collated,
                    logger=logger)

    fill = textwrap.fill(f'{"[INFO]:":10} Finished extracting putative ortholog trees using the Monophyletic '
                         f'Outgroups (MO) '
                         f'algorithm. Trees have been written to directory: "{output_folder}".',
                         width=90, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)

    logger.info(f'{fill}')
