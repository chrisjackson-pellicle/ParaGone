#!/usr/bin/env python

# Adapted from Yang and Smith (2014) by Chris Jackson chris.jackson@rbg.vic.gov.au
# https://github.com/chrisjackson-pellicle

"""
Trim tips that sticking out (> relative_cutoff and >10 times longer than sister)
Also trim any tips that are > absolute_cutoff
"""

import os
import textwrap
import glob

from paragone import phylo3
from paragone.tree_utils import *
from paragone import utils


# return the outlier tip, with abnormal high contrast and long branch
def check_contrast_outlier(node0,
                           node1,
                           above0_branch_length,
                           above1_branch_length,
                           relative_cutoff,
                           tree_name=None,
                           logger=None):
    """
    Return the outlier tip, with abnormal high contrast and long branch, if present, else return None

    :param phylo3.Node node0: first node for comparison
    :param phylo3.Node node1: second node for comparison
    :param float above0_branch_length: branch length of first node
    :param float above1_branch_length: branch length of second node
    :param float relative_cutoff: relative cutoff for removing tree tips
    :param str tree_name: name of the tree e.g. 6886.paralogs.aln.hmm.trimmed.fasta.treefile
    :param logging.Logger logger: a logger object
    :return phylo3.Node/None: return outlier node if present, else None
    """

    if node0.istip and above0_branch_length > relative_cutoff:
        # CJJ I DON'T UNDERSTAND THIS - WHY COMPARE TO ZERO BRANCH LENGTH: IS THIS IF KINKS HAVE BEEN INTRODUCED?
        if above1_branch_length == 0.0:
            return node0, f'Branch length {above0_branch_length} is above relative cut-off value of ' \
                          f'{relative_cutoff}, and sister branch length is {above1_branch_length}'

        elif above0_branch_length / above1_branch_length > 10:  # i.e. more than ten times longer than sister branch
            return node0, f'Branch length {above0_branch_length} is above relative cut-off value of ' \
                          f'{relative_cutoff}, and is >10 times longer than sister branch ({above1_branch_length})'

    if node1.istip and above1_branch_length > relative_cutoff:
        # CJJ I DON'T UNDERSTAND THIS - WHY COMPARE TO ZERO BRANCH LENGTH, IS THIS IF KINKS HAVE BEEN INTRODUCED?
        if above0_branch_length == 0.0:
            return node1, f'Branch length {above1_branch_length} is above relative cut-off value of ' \
                          f'{relative_cutoff}, and sister branch length is {above0_branch_length}'

        elif above1_branch_length / above0_branch_length > 10:
            return node1, f'Branch length {above1_branch_length} is above relative cut-off value of ' \
                          f'{relative_cutoff}, and is >10 times longer than sister branch ({above0_branch_length})'

    return None, None


def remove_a_tip(root,
                 tip_node,
                 tree_name=None,
                 logger=None):
    """
    Remove a given tip from a tree, and remove the 'kink' that this produces.

    :param phylo3.Node root: tree object parsed by newick3.parse
    :param phylo3.Node tip_node:
    :param str tree_name: name of the tree e.g. 6886.paralogs.aln.hmm.trimmed.fasta.treefile
    :param logging.Logger logger: a logger object
    :return phylo3.Node/None root/None: pruned tree if more than four tips remaining, else None
    """

    node = tip_node.prune()

    if len(root.leaves()) > 3:
        node, root = remove_kink(node, root)
        return root
    else:
        logger.warning(f'{"[WARNING]:":10} After removing tip {tip_node.label}, tree {tree_name} has fewer than '
                       f'four tips left.')
        return None


def trim(curroot,
         relative_cutoff,
         absolute_cutoff,
         tree_name=None,
         logger=None):
    """
    Removes tree tips longer than specified cutoffs

    :param phylo3.Node curroot: tree object parsed by newick3.parse
    :param float relative_cutoff: relative cutoff for removing tree tips
    :param float absolute_cutoff: absolute cutoff for removing tree tips
    :param str tree_name: name of the tree e.g. 6886.paralogs.aln.hmm.trimmed.fasta.treefile
    :param logging.Logger logger: a logger object
    :return:
    """

    if curroot.nchildren == 2:
        print('YEAH')
        temp, root = remove_kink(curroot, curroot)  # CJJ not used?

    going = True
    nodes_above_absolute_cutoff = defaultdict(list)
    nodes_above_relative_cutoff = defaultdict(list)

    while going and curroot and len(curroot.leaves()) > 3:
        going = False
        for node in curroot.iternodes(order=1):  # POSTORDER

            # Check if node is a tip, and remove it if branch length is greater than absolute cutoff:
            if node.nchildren == 0:  # at the tip
                node.data['len'] = node.length
                # print(f'node.length is; {node.length}')
                # print(f'absolute_cutoff is: {absolute_cutoff}')
                if node.length > absolute_cutoff:
                    logger.debug(f'Tip {node.label} is on a branch with length {node.length}. Absolute cutoff is '
                                 f'{absolute_cutoff}. This tip will be removed.')
                    nodes_above_absolute_cutoff[node.label].extend([node.length, absolute_cutoff])
                    # print(nodes_above_absolute_cutoff)
                    curroot = remove_a_tip(curroot,
                                           node,
                                           tree_name=tree_name,
                                           logger=logger)
                    going = True
                    break
            elif node.nchildren == 1:  # kink in tree
                print('KINK')
                remove_kink(node, curroot)
                going = True
                break
            elif node.nchildren == 2:  # normal bifurcating internal nodes
                child0_node, child1_node = node.children[0], node.children[1]
                above0_branch_length, above1_branch_length = child0_node.data['len'], child1_node.data['len']
                # CJJ internal branch length is adjusted below - why?:
                node.data['len'] = ((above0_branch_length + above1_branch_length) / 2.) + node.length  # stepwise average
                outlier, reason = check_contrast_outlier(child0_node,
                                                         child1_node,
                                                         above0_branch_length,
                                                         above1_branch_length,
                                                         relative_cutoff,
                                                         tree_name=tree_name,
                                                         logger=logger)
                # if outlier != None:
                if outlier:
                    nodes_above_relative_cutoff[outlier.label].extend([reason, 'birfurcating internal node'])
                    curroot = remove_a_tip(curroot,
                                           outlier,
                                           tree_name=tree_name,
                                           logger=logger)
                    going = True  # need to keep checking
                    break
            else:  # 3 or more branches from this node. Pair-wise comparison
                total_len = 0
                nchild = node.nchildren
                for child in node.children:
                    total_len += child.data['len']
                node.data['len'] = total_len / float(node.nchildren)
                keep_checking = True
                for index1 in range(nchild):  # do all the pairwise comparison
                    for index2 in range(nchild):
                        if index2 <= index1:
                            continue  # avoid repeatedly checking a pair
                        child1_node, child2_node = node.children[index1], node.children[index2]
                        above1_branch_length, above2_branch_length = child1_node.data['len'], child2_node.data['len']
                        outlier, reason = check_contrast_outlier(child1_node,
                                                                 child2_node,
                                                                 above1_branch_length,
                                                                 above2_branch_length,
                                                                 relative_cutoff,
                                                                 logger=logger)
                        if outlier:
                            nodes_above_relative_cutoff[outlier.label].extend([reason, '>=3 branches at node'])
                            curroot = remove_a_tip(curroot,
                                                   outlier,
                                                   tree_name=tree_name,
                                                   logger=logger)
                            going = True  # need to keep checking
                            keep_checking = False  # to break the nested loop
                            break
                    if not keep_checking:
                        break

    return curroot, nodes_above_absolute_cutoff, nodes_above_relative_cutoff


def write_trim_report(collated_trim_report_dict,
                      report_directory,
                      logger=None):
    """
    Writes a *.tsv report detailing which tips were trimmed from each tree, and why.

    :param dict collated_trim_report_dict: dictionary of default dicts for absolute and relative cut-off tips/reasons
    :param str report_directory: path to directory for report files
    :param logging.Logger logger: a logger object
    :return:
    """

    report_filename = f'{report_directory}/trees_trimmed_report.tsv'

    logger.info('')
    fill = utils.fill_forward_slash(f'{"[INFO]:":10} Writing trim tips report to file: "{report_filename}"',
                                    width=90, subsequent_indent=' ' * 11, break_on_forward_slash=True)

    logger.info(f'{fill}')

    all_tree_stats_for_report = []

    for tree_name, dictionaries in collated_trim_report_dict.items():

        tree_stats = [tree_name]

        try:
            check = dictionaries['absolute_cutoff']
            assert len(check) != 0
            tree_stats.append(len(check))  # This doesn't record actual lengths
            tips = [key for key in check.keys()]
            tree_stats.append('; '.join(tips))
        except AssertionError:
            tree_stats.append('0')
            tree_stats.append('N/A')

        try:
            check = dictionaries['relative_cutoff']
            assert len(check) != 0
            tree_stats.append(len(check))  # This doesn't record actual lengths
            tips = [key for key in check.keys()]
            tree_stats.append('; '.join(tips))
        except AssertionError:
            tree_stats.append('0')
            tree_stats.append('N/A')

        try:
            check = dictionaries['trimmed_trees_greater_than_four_taxa']
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries['trimmed_trees_fewer_than_four_taxa']
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        all_tree_stats_for_report.append(tree_stats)

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'Tree name\t'
                            f'Tips removed > absolute cutoff\t'
                            f'Tip names\t'
                            f'Tips removed > relative cutoff\t'
                            f'Tip names\t'
                            f'Trimmed trees > four taxa\t'
                            f'Trimmed trees < four taxa'
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

    logger.debug(f'{"[INFO]:":10} Module trim_tree_tips was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    logger.info(f'{"[INFO]:":10} ======> TRIMMING TREE TIPS <======\n')

    logger.debug(f'{"[INFO]:":10} Module trim_tree_tips was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    # Checking input directories and files:
    tree_file_directory = '05_trees_pre_quality_control'
    tree_file_suffix = '.treefile'

    directory_suffix_dict = {tree_file_directory: tree_file_suffix}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Create output folder for trimmed trees:
    treefile_directory_basename = os.path.basename(tree_file_directory)
    output_folder = f'06_{treefile_directory_basename.lstrip("05_")}_trimmed'
    trimmed_tree_output_folder = utils.createfolder(output_folder)

    collated_trim_report_dict = defaultdict(lambda: defaultdict())

    for treefile in glob.glob(f'{tree_file_directory}/*{tree_file_suffix}'):
        basename = os.path.basename(treefile)
        with open(treefile, 'r') as treefile_handle:
            intree = newick3.parse(treefile_handle.readline())

        logger.info(f'{"[INFO]:":10} Analysing tree: {basename}')

        trimmed_tree, nodes_above_absolute_cutoff, nodes_above_relative_cutoff = \
            trim(intree,
                 float(args.trim_tips_relative_cutoff),
                 float(args.trim_tips_absolute_cutoff),
                 tree_name=basename,
                 logger=logger)

        collated_trim_report_dict[basename]['absolute_cutoff'] = nodes_above_absolute_cutoff
        collated_trim_report_dict[basename]['relative_cutoff'] = nodes_above_relative_cutoff

        if trimmed_tree:
            collated_trim_report_dict[basename]['trimmed_trees_greater_than_four_taxa'] = trimmed_tree
            with open(f'{output_folder}/{basename}.tt', 'w') as outfile_handle:
                outfile_handle.write(newick3.tostring(trimmed_tree) + ";\n")
        else:
            collated_trim_report_dict[basename]['trimmed_trees_fewer_than_four_taxa'] = trimmed_tree
            logger.warning(f'{"[WARNING]:":10} No trimmed tree produced for {basename}!')

    # Write a report of tips trimmed from each tree, and why:
    write_trim_report(collated_trim_report_dict,
                      report_directory,
                      logger=logger)

    fill = textwrap.fill(f'{"[INFO]:":10} Finished trimming tips of input trees. Trimmed trees have been written to '
                         f'directory: "{output_folder}".',
                         width=90, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)

    logger.info(f'{fill}')

