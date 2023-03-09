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
import subprocess
import shutil
from ete3 import Tree

from paragone.tree_utils import *
from paragone import utils


def write_trim_report(collated_trim_report_dict,
                      report_directory,
                      min_tips,
                      quantile,
                      logger=None):
    """
    Writes a *.tsv report detailing which tips were trimmed from each tree, and why.

    :param dict collated_trim_report_dict: dictionary of default dicts for absolute and relative cut-off tips/reasons
    :param str report_directory: path to directory for report files
    :param int min_tips: the minimum number of tips in a tree after trimming tips
    :param float quantile: quantile value provided to TreeShrink
    :param logging.Logger logger: a logger object
    :return:
    """

    report_filename = f'{report_directory}/trees_trimmed_report.tsv'

    logger.info('')
    fill = utils.fill_forward_slash(f'{"[INFO]:":10} Writing trim tips report to file: "{report_filename}"',
                                    width=90, subsequent_indent=' ' * 11, break_on_forward_slash=True)

    logger.info(f'{fill}')

    all_tree_stats_for_report = []

    for tree_name, dictionaries in sorted(collated_trim_report_dict.items()):

        tree_stats = [tree_name]

        try:
            check = dictionaries['tips_removed']
            assert len(check) != 0
            tree_stats.append(len(check))  # This doesn't record actual lengths
            tips = [tip for tip in check]
            tree_stats.append('; '.join(tips))
        except AssertionError:
            tree_stats.append('0')
            tree_stats.append('N/A')

        try:
            check = dictionaries[f'trimmed_trees_greater_than_{min_tips}_taxa']
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries[f'trimmed_trees_fewer_than_{min_tips}_taxa']
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        all_tree_stats_for_report.append(tree_stats)

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'Tree name\t'
                            f'Tips removed by TreeShrink with quantile {quantile}\t'
                            f'Tip names\t'
                            f'Trimmed trees > {min_tips} taxa\t'
                            f'Trimmed trees < {min_tips} taxa'
                            f'\n')

        for stats in all_tree_stats_for_report:
            stats_joined = '\t'.join([str(stat) for stat in stats])
            report_handle.write(f'{stats_joined}\n')


def treeshrink(treefile,
               trimmed_tree_output_folder,
               q_value=0.05,
               logger=None):
    """
    Runs the program TreeShrink (https://github.com/uym2/TreeShrink) using the 'per-gene' model on a single tree file.

    :param str treefile: name of tree file in newick format
    :param str trimmed_tree_output_folder: name of output folder for trimmed trees
    :param float q_value: value provided the TreeShrink --quantile parameter; default is 0.05
    :param logging.Logger logger: a logger object
    :return str, list: expected_shrunk_treefile_renamed, tips_removed: name of (renamed and moved) tree file output
    by TreeShrink, and a list of tips removed (if any).
    """

    treeshrink_output_archive = utils.createfolder(f'{trimmed_tree_output_folder}/treeshrink_output_archive')

    treefile_basename = os.path.basename(treefile)
    expected_shrunk_treefile_dir = f'{trimmed_tree_output_folder}/{treefile_basename}.ts_dir'
    expected_shrunk_treefile = f'{trimmed_tree_output_folder}/{treefile_basename}.ts_dir/output.treefile'
    expected_shrunk_treefile_renamed = f'{trimmed_tree_output_folder}/{treefile_basename}.tt'

    logger.info(f'{"[INFO]:":10} Running TreeShrink for tree {treefile_basename}')

    try:
        assert utils.file_exists_and_not_empty(expected_shrunk_treefile_renamed)
        logger.debug(f'TreeShrink output file {expected_shrunk_treefile_renamed} exists for {treefile_basename}, '
                     f'skipping...')

        # Get the names of removed tips, if any:
        with open(f'{treeshrink_output_archive}/{treefile_basename}.ts_dir/output.txt', 'r') as output_handle:
            tips_removed = []
            for line in output_handle.readlines():
                if line.strip():
                    tip_names = line.strip().split('\t')
                    tips_removed.extend(tip_names)

        return expected_shrunk_treefile_renamed, tips_removed

    except AssertionError:

        try:
            shutil.rmtree(expected_shrunk_treefile_dir)  # --force option to TreeShrink appends trees to existing file
        except FileNotFoundError:
            pass

        treeshrink_command = f'run_treeshrink.py --tree {treefile} --centroid --force --mode per-gene --quantiles' \
                             f' {str(q_value)} --outdir {trimmed_tree_output_folder}/{treefile_basename}.ts_dir'

        logger.debug(f'Running command: {treeshrink_command}')

        try:
            result = subprocess.run(treeshrink_command,
                                    universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    shell=True, check=True)
            logger.debug(f'treeshrink check_returncode() is: {result.check_returncode()}')
            logger.debug(f'treeshrink stdout is: {result.stdout}')
            logger.debug(f'treeshrink stderr is: {result.stderr}')

            # Get the names of removed tips, if any:
            with open(f'{expected_shrunk_treefile_dir}/output.txt', 'r') as output_handle:
                tips_removed = []
                for line in output_handle.readlines():
                    if line.strip():
                        tip_names = line.strip().split('\t')
                        tips_removed.extend(tip_names)

        except subprocess.CalledProcessError as exc:
            logger.error(f'treeshrink FAILED. Output is: {exc}')
            logger.error(f'v stdout is: {exc.stdout}')
            logger.error(f'treeshrink stderr is: {exc.stderr}')
            raise ValueError('There was an issue running TreeShrink. Check input files!')

        # Unroot the TreeShrink tree and write it to the parent directory:
        try:

            logger.debug(f'Writing shrunk tree file to file {expected_shrunk_treefile_renamed}')
            tree = Tree(expected_shrunk_treefile)
            tree.unroot()
            tree.write(outfile=f'{expected_shrunk_treefile_renamed}', format=0)

            # Move TreeShrink directory to folder treeshrink_output_dirs
            shutil.move(expected_shrunk_treefile_dir, treeshrink_output_archive)

        except FileNotFoundError:
            raise

        return expected_shrunk_treefile_renamed, tips_removed


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

    logger.debug(f'{"[INFO]:":10} Module trim_tree_treeshrink was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    logger.info(f'\n{"[INFO]:":10} ======> TRIMMING TREE TIPS USING TREESHRINK<======\n')

    # Checking input directories and files:
    tree_file_directory = '05_trees_pre_quality_control'
    tree_file_suffix = '.treefile'

    directory_suffix_dict = {tree_file_directory: tree_file_suffix}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Create output folder for trimmed trees:
    output_folder = f'06_trees_trimmed'
    trimmed_tree_output_folder = utils.createfolder(output_folder)

    collated_trim_report_dict = defaultdict(lambda: defaultdict())

    for treefile in sorted(glob.glob(f'{tree_file_directory}/*{tree_file_suffix}')):
        basename = os.path.basename(treefile)
        trimmed_tree, tips_removed = treeshrink(treefile,
                                                trimmed_tree_output_folder,
                                                q_value=args.treeshrink_q_value,
                                                logger=logger)

        collated_trim_report_dict[basename]['tips_removed'] = tips_removed

        # Read in the tree produced by TreeShrink:
        with open(trimmed_tree, 'r') as trimmed_tree_handle:
            intree = newick3.parse(trimmed_tree_handle.readline())

        # Check number of tips:
        if len(intree.leaves()) >= args.min_tips:
            collated_trim_report_dict[basename][f'trimmed_trees_greater_than_{args.min_tips}_taxa'] = intree

        else:
            collated_trim_report_dict[basename][f'trimmed_trees_fewer_than_{args.min_tips}_taxa'] = intree
            logger.warning(f'{"[WARNING]:":10} No trimmed tree produced for {basename}!')
            os.remove(trimmed_tree)

    # Write a report of tips trimmed from each tree, and why:
    write_trim_report(collated_trim_report_dict,
                      report_directory,
                      args.min_tips,
                      args.treeshrink_q_value,
                      logger=logger)

    fill = textwrap.fill(f'{"[INFO]:":10} Finished trimming tips of input trees. Trimmed trees have been written to '
                         f'directory: "{output_folder}".',
                         width=90, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)

    logger.info(f'{fill}')
