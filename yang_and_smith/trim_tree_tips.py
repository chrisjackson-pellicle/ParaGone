#!/usr/bin/env python

# Adapted from Yang and Smith by Chris Jackson chris.jackson@rbg.vic.gov.au

"""
Trim tips that sticking out (> relative_cutoff and >10 times longer than sister)
Also trim any tips that are > absolute_cutoff
"""

import os
import sys
import textwrap
import glob
from collections import defaultdict
import newick3
import phylo3
from tree_utils import *

from yang_and_smith import utils


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
        # CJJ I DON'T UNDERSTAND THIS - WHY COMPARE TO ZERO BRANCH LENGTH, IS THIS IF KINKS HAVE BEEN INTRODUCED?
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


def remove_a_tip(root, tip_node, tree_name=None, logger=None):
    """
    Remove a given tip from a tree, and remove the 'kink' that this produces.

    :param phylo3.Node root:
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
        logger.warning(f'{"[WARNING]:":10} After removing tip {tip_node.label}, tree {tree_name} has less than four '
                       f'tips left.')
        return None


def trim(curroot, relative_cutoff, absolute_cutoff, tree_name=None, logger=None):
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
                if node.length > absolute_cutoff:
                    logger.debug(f'Tip {node.label} is on a branch with length {node.length}. Absolute cutoff is '
                                 f'{absolute_cutoff}. This tip will be removed.')
                    nodes_above_absolute_cutoff[node.label].extend([node.length, absolute_cutoff])
                    curroot = remove_a_tip(curroot, node)
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


def write_trim_report(collated_trim_report_dict, treefile_directory, logger=None):
    """
    Writes a *.tsv report detailing which tips were trimmed from each tree, and why.

    :param dict collated_trim_report_dict: dictionary of default dicts for absolute and relative cut-off tips/reasons
    :param str treefile_directory: name of tree file directory for report filename
    :param logging.Logger logger: a logger object
    :return:
    """

    basename = os.path.basename(treefile_directory)
    report_filename = f'{basename}_trim_tips_report.tsv'

    logger.info(f'{"[INFO]:":10} Writing trim tips report to file {report_filename}')

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'Tree_name\tTips_removed\tTip_name\tReason_for_removal\tFrom_node_type\n')
        for tree_name, dictionaries in collated_trim_report_dict.items():
            absolute_dict = dictionaries['absolute_cutoff']
            relative_dict = dictionaries['relative_cutoff']

            if not absolute_dict and not relative_dict:
                report_handle.write(f'{tree_name}\tN\tN/A\tN/A\tN/A\n')

            for trimmed_tip, data in absolute_dict.items():
                branch_length, absolute_cutoff = data
                report_handle.write(f'{tree_name}\tY\t{trimmed_tip}\tBranch length {branch_length} above absolute '
                                    f'cut-off value {absolute_cutoff}\tN/A\n')

            for trimmed_tip, data in relative_dict.items():
                reason, node_type = data
                report_handle.write(f'{tree_name}\tY\t{trimmed_tip}\t{reason}\t{node_type}\n')


def main(args):
    """
    Entry point for the resolve_paralogs.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = utils.setup_logger(__name__, 'trim_tree_tips')

    logger.info(f'{"[INFO]:":10} Subcommand trim_tree_tips was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')
    logger.info(f'{"[INFO]:":10} Relative cutoff value: {args.relative_cutoff}')
    logger.info(f'{"[INFO]:":10} Absolute cutoff value: {args.absolute_cutoff}')

    filecount = 0
    utils.createfolder(args.output_folder)

    collated_trim_report_dict = dict()

    for treefile in glob.glob(f'{args.treefile_directory}/*{args.tree_file_suffix}'):
        basename = os.path.basename(treefile)
        filecount += 1
        with open(treefile, 'r') as treefile_handle:
            intree = newick3.parse(treefile_handle.readline())

        logger.info(f'{"[INFO]:":10} Analysing tree: {basename}')

        trimmed_tree, nodes_above_absolute_cutoff, nodes_above_relative_cutoff = \
            trim(intree,
                 float(args.relative_cutoff),
                 float(args.absolute_cutoff),
                 tree_name=basename,
                 logger=logger)

        collated_trim_report_dict[basename] = {'absolute_cutoff': nodes_above_absolute_cutoff,
                                               'relative_cutoff': nodes_above_relative_cutoff}

        if trimmed_tree:
            with open(f'{args.output_folder}/{basename}.tt', 'w') as outfile_handle:
                outfile_handle.write(newick3.tostring(trimmed_tree) + ";\n")
        else:
            logger.warning(f'No trimmed tree produced for {basename}!')

    # Write a report of tips trimmed from each tree, and why:
    write_trim_report(collated_trim_report_dict,
                      args.treefile_directory,
                      logger=logger)

    assert filecount > 0, f'No files with suffix {args.tree_file_suffix} found in {args.treefile_directory}'
