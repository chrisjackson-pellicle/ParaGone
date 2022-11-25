#!/usr/bin/env python

# Adapted from Yang and Smith (2014) by Chris Jackson chris.jackson@rbg.vic.gov.au
# https://github.com/chrisjackson-pellicle

"""
Input is a directory of trees that end with ".tt". This can be changed.

- Mask both mono- and paraphyletic (the latter optional) tips that belong to the same taxon
- Keeps the tip that has the greatest number of unambiguous characters in the QC'd alignment
"""

import os
import sys
import textwrap
import glob
from collections import defaultdict

from paragone import newick3
from paragone import phylo3
from paragone.tree_utils import get_name, remove_kink
from paragone.seq import read_fasta_file
from paragone import utils


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
    else None, default dictionary (list) of removed tips and associated data
    """

    going = True

    pruned_monophyletic_tips_dict = defaultdict(list)
    trees_with_fewer_than_four_tips_dict = {}

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
                        trees_with_fewer_than_four_tips_dict[tree_name] = curroot
                        return None, pruned_monophyletic_tips_dict, trees_with_fewer_than_four_tips_dict

                    going = True
                    break

    return curroot, pruned_monophyletic_tips_dict, trees_with_fewer_than_four_tips_dict


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
    trees_with_fewer_than_four_tips_dict = {}

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
                        trees_with_fewer_than_four_tips_dict[tree_name] = curroot
                        return None, pruned_paraphyletic_tips_dict, trees_with_fewer_than_four_tips_dict

                    going = True
                    break

    return curroot, pruned_paraphyletic_tips_dict


def write_mask_report(collated_mask_report_dict,
                      treefile_directory,
                      logger=None):
    """
    Writes a *.tsv report detailing which tips were masked (removed) from each tree, and why.

    :param dict collated_mask_report_dict: dictionary of default dicts for mono and (optional) paraphyletic cut-off
    tips/data
    :param str treefile_directory: name of tree file directory for report filename
    :param logging.Logger logger: a logger object
    :return:
    """

    basename = os.path.basename(treefile_directory)
    report_filename = f'00_logs_and_reports_resolve_paralogs/reports/{basename.lstrip("08_")}_masked_report.tsv'

    logger.info(f'{"[INFO]:":10} Writing mask tips report to file {report_filename}')

    all_tree_stats_for_report = []

    for tree_name, dictionaries in collated_mask_report_dict.items():

        tree_stats = [tree_name]

        try:
            check = dictionaries['monophyletic_tips']
            assert len(check) != 0
            tree_stats.append(len(check))  # This doesn't record actual unambiguous character counts
            tips = [key for key in check.keys()]
            tree_stats.append('; '.join(tips))
        except AssertionError:
            tree_stats.append('0')
            tree_stats.append('N/A')

        try:
            check = dictionaries['paraphyletic_tips']
            assert len(check) != 0
            tree_stats.append(len(check))  # This doesn't record actual unambiguous character counts
            tips = [key for key in check.keys()]
            tree_stats.append('; '.join(tips))
        except KeyError:
            tree_stats.append('0')
            tree_stats.append('N/A')
        except AssertionError:
            tree_stats.append('0')
            tree_stats.append('N/A')

        try:
            check = dictionaries['less_than_four_taxa_after_masking_mono']
            assert len(check) != 0
            tree_stats.append('Y')
        except AssertionError:
            tree_stats.append('N')

        try:
            check = dictionaries['less_than_four_taxa_after_masking_para']
            assert len(check) != 0
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')
        except AssertionError:
            tree_stats.append('N')

        all_tree_stats_for_report.append(tree_stats)

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'Tree name\t'
                            f'Monophyletic tips removed ("masked")\t'
                            f'Retained tip names\t'
                            f'Paraphyletic tips removed ("masked")\t'
                            f'Retained tip names\t'
                            f'Masked trees < four taxa after masking mono\t'
                            f'Masked trees < four taxa after masking para'
                            f'\n')

        for stats in all_tree_stats_for_report:
            stats_joined = '\t'.join([str(stat) for stat in stats])
            report_handle.write(f'{stats_joined}\n')


def main(args):
    """
    Entry point for the paragone_main.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = utils.setup_logger(__name__, '00_logs_and_reports_resolve_paralogs/logs/06_mask_tree_tips')

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
    logger.debug(args)

    # Checking input directories and files:
    directory_suffix_dict = {args.treefile_directory: args.tree_file_suffix,
                             args.alignment_directory: args.alignment_file_suffix}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Create output folder:
    treefile_directory_basename = os.path.basename(args.treefile_directory)
    output_folder = f'09_{treefile_directory_basename.lstrip("08_")}_masked'
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

    collated_mask_report_dict = defaultdict(lambda: defaultdict())

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
        curroot, \
        pruned_monophyletic_tips_dict, \
        trees_with_fewer_than_four_tips_mono_dict = \
            mask_monophyletic_tips(intree,
                                   unamiguous_characters_dict,
                                   tree_name=tree_file_basename,
                                   logger=logger)

        collated_mask_report_dict[tree_file_basename]['monophyletic_tips'] = pruned_monophyletic_tips_dict
        collated_mask_report_dict[tree_file_basename]['less_than_four_taxa_after_masking_mono'] = \
            trees_with_fewer_than_four_tips_mono_dict

        if not curroot:
            logger.warning(f'No masked tree produced for {tree_file_basename}!')
            continue

        # Optionally, remove ('mask' in Yang and Smith terminology) all but one tip for paraphyletic clades
        # from same taxon name:
        if args.remove_paraphyletic_tips:
            curroot, \
            pruned_paraphyletic_tips_dict, trees_with_fewer_than_four_tips_para_dict\
                = mask_paraphyletic_tips(curroot,
                                         unamiguous_characters_dict,
                                         tree_name=tree_file_basename,
                                         logger=logger)

            collated_mask_report_dict[tree_file_basename]['paraphyletic_tips'] = pruned_paraphyletic_tips_dict
            collated_mask_report_dict[tree_file_basename]['less_than_four_taxa_after_masking_para'] = \
                trees_with_fewer_than_four_tips_para_dict

        if not curroot:
            logger.warning(f'No masked tree produced for {tree_file_basename}!')
            continue

        # Write 'masked' tree to file:
        tree_to_write = f'{output_folder}/{tree_file_basename}.mm'
        with open(tree_to_write, "w") as tree_outfile:
            tree_outfile.write(newick3.tostring(curroot) + ";\n")

    fill = textwrap.fill(f'{"[INFO]:":10} Finished masking tips of input trees. Masked trees have been written to '
                         f'directory: "{output_folder}".',
                         width=90, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)

    logger.info(f'{fill}')

    # Write a report of tips trimmed from each tree, and why:
    write_mask_report(collated_mask_report_dict,
                      args.treefile_directory,
                      logger=logger)

    try:
        assert filecount > 0
    except AssertionError:
        logger.error(f'{"[ERROR]:":10} No files with suffix {args.tree_file_suffix} found in'
                     f' {args.treefile_directory}.')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Finished masking tree tips.')

