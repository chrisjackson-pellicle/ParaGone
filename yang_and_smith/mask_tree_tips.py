#!/usr/bin/env python

# Adapted from Yang and Smith by Chris Jackson chris.jackson@rbg.vic.gov.au

"""
Input is a directory of trees that end with ".tt". This can be changed.

- Mask both mono- and paraphyletic (optional) tips that belong to the same taxon
- Keeps the tip that has the most unambiguous characters in the QC'd alignment
"""

import os
import sys
import textwrap
import glob
from collections import defaultdict
import newick3
import phylo3
from tree_utils import get_name, remove_kink
from seq import read_fasta_file

from yang_and_smith import utils


def mask_monophyletic_tips(curroot,
                           unamiguous_characters_dict,
                           tree_name=None,
                           logger=None):
    """
    Removes all but one tip for a given taxon if monophyletic with same taxon name. Keep tip for which the
    corresponding sequence has the most unambiguous characters.

    :param phylo3.Node curroot: tree object parsed by newick3.parse
    :param dict unamiguous_characters_dict: dictionary of sequence name to number of unambiguous characters
    :param str tree_name: file name of the tree being examined
    :param logging.Logger logger: a logger object
    :return phylo3.Node/None,collections.defaultdict  curroot/None,pruned_monophyletic_tips_dict : pruned tree object of more than four tips,
    else None, default dictionary (list) of removed tips and associate data
    """

    going = True

    pruned_monophyletic_tips_dict = defaultdict(list)

    while going and len(curroot.leaves()) >= 4:
        going = False
        for node in curroot.iternodes():  # walk through nodes
            if not node.istip:
                continue  # only look at tips

            for sister in node.get_sisters():
                if sister.istip and get_name(node.label) == get_name(sister.label):  # masking
                    if unamiguous_characters_dict[node.label] > unamiguous_characters_dict[sister.label]:

                        pruned_monophyletic_tips_dict[sister.label].extend(
                            [node.label,
                             unamiguous_characters_dict[node.label],
                             unamiguous_characters_dict[sister.label]])

                        node = sister.prune()

                    else:
                        pruned_monophyletic_tips_dict[node.label].extend(
                            [sister.label,
                             unamiguous_characters_dict[sister.label],
                             unamiguous_characters_dict[node.label]])

                        node = node.prune()

                    if len(curroot.leaves()) >= 4:
                        if (node == curroot and node.nchildren == 2) or (node != curroot and node.nchildren == 1):
                            node, curroot = remove_kink(node, curroot)
                    else:
                        logger.warning(f'{"[WARNING]:":10} Tree {tree_name} has fewer than four tips after removal of '
                                       f'monophyletic tips!')
                        return None, pruned_monophyletic_tips_dict

                    going = True
                    break

    return curroot, pruned_monophyletic_tips_dict


def mask_paraphyletic_tips(curroot,
                           unamiguous_characters_dict,
                           tree_name=None,
                           logger=None):
    """
    Removes all but one tip for a given taxon if paraphyletic with same taxon name. Keep tip for which the
    corresponding sequence has the most unambiguous characters.

    :param phylo3.Node curroot: tree object parsed by newick3.parse
    :param dict unamiguous_characters_dict: dictionary of sequence name to number of unambiguous characters
    :param str tree_name: file name of the tree being examined
    :param logging.Logger logger: a logger object
    :return:
    """

    going = True

    pruned_paraphyletic_tips_dict = defaultdict(list)

    while going and len(curroot.leaves()) >= 4:
        going = False
        for node in curroot.iternodes():  # walk through nodes
            if not node.istip:
                continue  # only look at tips

            parent = node.parent
            if node == curroot or parent == curroot:
                continue  # no paraphyletic tips for the root

            for para in parent.get_sisters():
                # print(f'sister is: {newick3.tostring(para)}')
                if para.istip and get_name(node.label) == get_name(para.label):
                    if unamiguous_characters_dict[node.label] > unamiguous_characters_dict[para.label]:

                        pruned_paraphyletic_tips_dict[para.label].extend(
                            [node.label,
                             unamiguous_characters_dict[node.label],
                             unamiguous_characters_dict[para.label]])

                        node = para.prune()

                    else:
                        pruned_paraphyletic_tips_dict[node.label].extend(
                            [para.label,
                             unamiguous_characters_dict[para.label],
                             unamiguous_characters_dict[node.label]])

                        node = node.prune()

                    if len(curroot.leaves()) >= 4:
                        if (node == curroot and node.nchildren == 2) or (node != curroot and node.nchildren == 1):
                            node, curroot = remove_kink(node, curroot)
                    else:
                        logger.warning(f'{"[WARNING]:":10} Tree {tree_name} has fewer than four tips after removal of '
                                       f'paraphyletic tips!')
                        return None, pruned_paraphyletic_tips_dict

                    going = True
                    break

    return curroot, pruned_paraphyletic_tips_dict


def write_mask_report(collated_mask_report_dict, treefile_directory, logger=None):
    """
    Writes a *.tsv report detailing which tips were masked (removed) from each tree, and why.

    :param dict collated_mask_report_dict: dictionary of default dicts for mono and (optional) paraphyletic cut-off
    tips/data
    :param str treefile_directory: name of tree file directory for report filename
    :param logging.Logger logger: a logger object
    :return:
    """

    basename = os.path.basename(treefile_directory)
    report_filename = f'{basename}_masked_report.tsv'

    logger.info(f'{"[INFO]:":10} Writing mask tips report to file {report_filename}')

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'Tree_name\tTips_removed\tTip_name\tReason_for_removal\tFrom_node_type\n')
        for tree_name, dictionaries in collated_mask_report_dict.items():
            monophyletic_dict = dictionaries['monophyletic_tips']
            try:
                paraphyletic_dict = dictionaries['paraphyletic_tips']
            except KeyError:
                logger.debug(f'No paraphyletic tips dictionary found in {dictionaries}')
                paraphyletic_dict = None

            if not monophyletic_dict and not paraphyletic_dict:
                report_handle.write(f'{tree_name}\tN\tN/A\tN/A\tN/A\n')

            for trimmed_tip, data in monophyletic_dict.items():
                monopyletic_sister_retained,\
                monopyletic_sister_retained_character_count,\
                trimmed_tip_character_count \
                    = data

                report_handle.write(f'{tree_name}\tY\t{trimmed_tip}\tMonophyletic sister '
                                    f'{monopyletic_sister_retained} unambiguous character count is'
                                    f' {monopyletic_sister_retained_character_count}, whereas {trimmed_tip} '
                                    f'unambiguous character count is {trimmed_tip_character_count}\tMonophyletic\n')

            if paraphyletic_dict:
                for trimmed_tip, data in paraphyletic_dict.items():
                    parapyletic_sister_retained, \
                    parapyletic_sister_retained_character_count, \
                    trimmed_tip_character_count \
                        = data

                    report_handle.write(f'{tree_name}\tY\t{trimmed_tip}\tParaphyletic sister '
                                        f'{parapyletic_sister_retained} unambiguous character count is'
                                        f' {parapyletic_sister_retained_character_count}, whereas {trimmed_tip} '
                                        f'unambiguous character count is {trimmed_tip_character_count}\tParaphyletic\n')


def main(args):
    """
    Entry point for the resolve_paralogs.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = utils.setup_logger(__name__, 'logs_resolve_paralogs/05_mask_tree_tips')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand mask_tree_tips was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')

    # Create output folder:
    treefile_directory_basename = os.path.basename(args.treefile_directory)
    output_folder = f'{treefile_directory_basename}_masked'
    utils.createfolder(output_folder)

    filecount = 0

    # Create a dictionary of gene name to alignment file name:
    gene_name_to_alignment_dict = {}  # key is gene name, value is the alignment file name
    for alignment_file in glob.glob(f'{args.alignment_directory}/*{args.alignment_file_suffix}'):
        alignment_basename = os.path.basename(alignment_file)
        alignment_gene_name = os.path.basename(alignment_file).split('.')[0]
        try:
            gene_name_to_alignment_dict[alignment_gene_name] = alignment_basename
        except KeyError:
            sys.exit(f'Gene name {alignment_gene_name} occurs more than once in the directory'
                     f' {args.alignment_directory}')
        except:
            raise

    collated_mask_report_dict = dict()

    for tree_file in glob.glob(f'{args.treefile_directory}/*{args.tree_file_suffix}'):
        tree_file_basename = os.path.basename(tree_file)
        with open(tree_file, 'r') as tree_file_handle:
            intree = newick3.parse(tree_file_handle.readline())
        tree_gene_name = os.path.basename(tree_file).split('.')[0]
        filecount += 1

        logger.info(f'{"[INFO]:":10} Analysing tree: {tree_file_basename}')

        unamiguous_characters_dict = {}  # key is seqid, value is number of unambiguous chrs
        for seq in read_fasta_file(f'{args.alignment_directory}/{gene_name_to_alignment_dict[tree_gene_name]}'):
            seq.name = seq.name.split()[0]
            # The line above is necessary to remove extraneous information that otherwise gets captured in the
            # dictionary key below. This wouldn't happen with Biopython...
            seq.description = ''
            for ambiguous_character in ['-', 'X', "x", "?", "*"]:
                seq.seq = seq.seq.replace(ambiguous_character, '')  # ignore gaps, xs, Xs question marks and stop codons
            unamiguous_characters_dict[seq.name] = len(seq.seq)

        # Remove ('mask' in Yang and Smith terminology) all but one tip for monophyletic clades from same taxon name.
        # This is designed to filter alleles or close paralogs:
        curroot, pruned_monophyletic_tips_dict = mask_monophyletic_tips(intree,
                                                                        unamiguous_characters_dict,
                                                                        tree_name=tree_file_basename,
                                                                        logger=logger)

        collated_mask_report_dict[tree_file_basename] = {'monophyletic_tips': pruned_monophyletic_tips_dict}

        if not curroot:
            logger.warning(f'No masked tree produced for {tree_file_basename}!')
            continue

        # Optionally, remove ('mask' in Yang and Smith terminology) all but one tip for paraphyletic clades
        # from same taxon name:
        if args.remove_paraphyletic_tips:
            curroot, pruned_paraphyletic_tips_dict = mask_paraphyletic_tips(curroot,
                                                                            unamiguous_characters_dict,
                                                                            tree_name=tree_file_basename,
                                                                            logger=logger)

            collated_mask_report_dict[tree_file_basename]['paraphyletic_tips'] = pruned_paraphyletic_tips_dict

        if not curroot:
            logger.warning(f'No masked tree produced for {tree_file_basename}!')
            continue

        # Write 'masked' tree to file:
        tree_to_write = f'{output_folder}/{tree_file_basename}.mm'
        with open(tree_to_write, "w") as tree_outfile:
            tree_outfile.write(newick3.tostring(curroot) + ";\n")

    # Write a report of tips trimmed from each tree, and why:
    write_mask_report(collated_mask_report_dict,
                      args.treefile_directory,
                      logger=logger)

    assert filecount > 0, logger.error(f'{"[ERROR]:":10} No files with suffix {args.tree_file_suffix} found in'
                                       f' {args.treefile_directory}')


