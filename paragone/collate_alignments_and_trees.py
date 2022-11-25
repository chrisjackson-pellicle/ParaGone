#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
- Collates all HmmCleaned alignments (i.e. in folders starting with prefix '04_batch' in to a single folder.
- Collate all newick tree files (i.e. in folders starting with prefix '05_batch' in to a single folder.
"""
import shutil
import textwrap
import glob
import sys
import os
import re

from paragone import utils


def collate_alignment_to_tree(directory_contents,
                              logger=None):
    """
    Collates all alignments and tree from the step "alignment_to_tree".

    :param list directory_contents: list of current working directory contents
    :param logging.Logger logger: a logger object
    :return:
    """

    output_alignment_directory = f'06_all_hmmcleaned_alignments'
    utils.createfolder(output_alignment_directory)

    output_treefile_directory = f'07_all_hmmcleaned_alignment_trees'
    utils.createfolder(output_treefile_directory)

    alignment_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
                             directory.startswith('04_batch')]

    treefile_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
                            directory.startswith('05_batch')]

    if not alignment_directories:
        logger.error(f'{"[ERROR]:":10} No alignment directories with names beginning "04_batch" found!')
        sys.exit(1)

    if not treefile_directories:
        logger.error(f'{"[ERROR]:":10} No tree file directories with names beginning "05_batch" found!')
        sys.exit(1)

    # Checking input directories and files:
    directory_suffix_dict = {}
    file_list = []
    for alignment_directory in alignment_directories:
        directory_suffix_dict[alignment_directory] = '.aln.hmm.trimmed.fasta'

    for treefile_directory in treefile_directories:
        directory_suffix_dict[treefile_directory] = '.treefile'

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Copy files to new directories:
    for alignment_directory in alignment_directories:
        for alignment_file in glob.glob(f'{alignment_directory}/*.aln.hmm.trimmed.fasta'):
            shutil.copy(alignment_file, output_alignment_directory)

    logger.info(f'{"[INFO]:":10} All HmmCleaned alignments have been copied to the directory: '
                f'"{output_alignment_directory}"')

    for treefile_directory in treefile_directories:
        for treefile in glob.glob(f'{treefile_directory}/*.treefile'):
            shutil.copy(treefile, output_treefile_directory)

    logger.info(f'{"[INFO]:":10} All trees from HmmCleaned alignments have been copied to the directory: '
                f'"{output_treefile_directory}"')


def collate_align_selected_and_tree(directory_contents,
                                    logger=None):
    """
    Collates all alignments and tree from the step "align_selected_and_tree".

    :param list directory_contents: list of current working directory contents
    :param logging.Logger logger: a logger object
    :return:
    """

    output_alignment_directory = f'16_all_trimmed_masked_cut_selected_alignments'
    utils.createfolder(output_alignment_directory)

    output_treefile_directory = f'17_all_trimmed_masked_cut_selected_alignments_trees'
    utils.createfolder(output_treefile_directory)

    alignment_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
                             directory.startswith('14_selected_batch')]

    treefile_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
                            directory.startswith('15_selected_batch')]

    if not alignment_directories:
        logger.error(f'{"[ERROR]:":10} No alignment directories with names beginning "04_batch" found!')
        sys.exit(1)

    if not treefile_directories:
        logger.error(f'{"[ERROR]:":10} No tree file directories with names beginning "05_batch" found!')
        sys.exit(1)

    # Checking input directories and files:
    directory_suffix_dict = {}
    file_list = []
    for alignment_directory in alignment_directories:
        directory_suffix_dict[alignment_directory] = '.aln.trimmed.fasta'

    for treefile_directory in treefile_directories:
        directory_suffix_dict[treefile_directory] = '.treefile'

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Copy files to new directories:
    for alignment_directory in alignment_directories:
        for alignment_file in glob.glob(f'{alignment_directory}/*.aln.trimmed.fasta'):
            shutil.copy(alignment_file, output_alignment_directory)

    logger.info(f'{"[INFO]:":10} All trimmed/masked/cut selected alignments have been copied to the directory: '
                f'"{output_alignment_directory}"')

    for treefile_directory in treefile_directories:
        for treefile in glob.glob(f'{treefile_directory}/*.treefile'):
            shutil.copy(treefile, output_treefile_directory)

    logger.info(f'{"[INFO]:":10} All trees from trimmed/masked/cut selected alignments have been copied to the '
                f'directory: "{output_treefile_directory}"')


def collate_pruned_alignments(algorithm_suffix,
                              directory_contents,
                              logger=None):
    """
    Collates all final alignments (names stripped, ready for concatentation) from given pruning algorithm.

    :param bool/str algorithm_suffix: if not False, name of pruning method mo/rt/mi
    :param list directory_contents: list of current working directory contents
    :param logging.Logger logger: a logger object
    :return:
    """

    output_alignment_directory = f'25_final_alignments_{algorithm_suffix}'
    utils.createfolder(output_alignment_directory)

    regex_string = f'.*{algorithm_suffix}_[0-9]+_stripped_names_alignments'
    regex = re.compile(regex_string)

    alignment_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
                             re.search(regex, str(directory))]

    if not alignment_directories:
        logger.error(f'{"[ERROR]:":10} No alignment directories matching pattern {regex_string} found!')
        sys.exit(1)

    # Checking input directories and files:
    directory_suffix_dict = {}
    file_list = []
    for alignment_directory in alignment_directories:
        directory_suffix_dict[alignment_directory] = 'stripped.aln.trimmed.fasta'

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Copy files to new directories:
    for alignment_directory in alignment_directories:
        for alignment_file in glob.glob(f'{alignment_directory}/*stripped.aln.trimmed.fasta'):
            shutil.copy(alignment_file, output_alignment_directory)

    logger.info(f'{"[INFO]:":10} All final alignments for algorithm {algorithm_suffix }have been copied to the '
                f'directory: "{output_alignment_directory}"')


def main(args):
    """
    Entry point for the paragone_main.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    if args.from_alignment_to_tree:
        logger = utils.setup_logger(__name__,
                                    '00_logs_and_reports_resolve_paralogs/logs/04_collate_alignments_and_trees')
    elif args.from_align_selected_and_tree:
        logger = utils.setup_logger(__name__,
                                    '00_logs_and_reports_resolve_paralogs/logs/10_collate_alignments_and_trees')

    elif args.from_prune_paralogs:
        logger = utils.setup_logger(__name__,
                                    f'00_logs_and_reports_resolve_paralogs/logs/'
                                    f'16_collate_final_alignments_{args.from_prune_paralogs}')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand collate_alignments_and_trees was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')
    logger.debug(args)

    directory_contents = os.listdir('.')

    # Create output folders:
    if args.from_alignment_to_tree:
        collate_alignment_to_tree(directory_contents,
                                  logger=logger)

        # output_alignment_directory = f'06_all_hmmcleaned_alignments'
        # utils.createfolder(output_alignment_directory)
        #
        # output_treefile_directory = f'07_all_hmmcleaned_alignment_trees'
        # utils.createfolder(output_treefile_directory)
        #
        # alignment_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
        #                          directory.startswith('04_batch')]
        #
        # treefile_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
        #                         directory.startswith('05_batch')]
        #
        # if not alignment_directories:
        #     logger.error(f'{"[ERROR]:":10} No alignment directories with names beginning "04_batch" found!')
        #     sys.exit(1)
        #
        # if not treefile_directories:
        #     logger.error(f'{"[ERROR]:":10} No tree file directories with names beginning "05_batch" found!')
        #     sys.exit(1)
        #
        # # Checking input directories and files:
        # directory_suffix_dict = {}
        # file_list = []
        # for alignment_directory in alignment_directories:
        #     directory_suffix_dict[alignment_directory] = '.aln.hmm.trimmed.fasta'
        #
        # for treefile_directory in treefile_directories:
        #     directory_suffix_dict[treefile_directory] = '.treefile'
        #
        # utils.check_inputs(directory_suffix_dict,
        #                    file_list,
        #                    logger=logger)
        #
        # # Copy files to new directories:
        # for alignment_directory in alignment_directories:
        #     for alignment_file in glob.glob(f'{alignment_directory}/*.aln.hmm.trimmed.fasta'):
        #         shutil.copy(alignment_file, output_alignment_directory)
        #
        # logger.info(f'{"[INFO]:":10} All HmmCleaned alignments have been copied to the directory: '
        #             f'"{output_alignment_directory}"')
        #
        # for treefile_directory in treefile_directories:
        #     for treefile in glob.glob(f'{treefile_directory}/*.treefile'):
        #         shutil.copy(treefile, output_treefile_directory)
        #
        # logger.info(f'{"[INFO]:":10} All trees from HmmCleaned alignments have been copied to the directory: '
        #             f'"{output_treefile_directory}"')

    elif args.from_align_selected_and_tree:
        collate_align_selected_and_tree(directory_contents,
                                        logger=logger)


    #     output_alignment_directory = f'16_all_trimmed_masked_cut_selected_alignments'
    #     utils.createfolder(output_alignment_directory)
    #
    #     output_treefile_directory = f'17_all_trimmed_masked_cut_selected_alignments_trees'
    #     utils.createfolder(output_treefile_directory)
    #
    #     alignment_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
    #                              directory.startswith('14_selected_batch')]
    #
    #     treefile_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
    #                             directory.startswith('15_selected_batch')]
    #
    #     if not alignment_directories:
    #         logger.error(f'{"[ERROR]:":10} No alignment directories with names beginning "04_batch" found!')
    #         sys.exit(1)
    #
    #     if not treefile_directories:
    #         logger.error(f'{"[ERROR]:":10} No tree file directories with names beginning "05_batch" found!')
    #         sys.exit(1)
    #
    #     # Checking input directories and files:
    #     directory_suffix_dict = {}
    #     file_list = []
    #     for alignment_directory in alignment_directories:
    #         directory_suffix_dict[alignment_directory] = '.aln.trimmed.fasta'
    #
    #     for treefile_directory in treefile_directories:
    #         directory_suffix_dict[treefile_directory] = '.treefile'
    #
    #     utils.check_inputs(directory_suffix_dict,
    #                        file_list,
    #                        logger=logger)
    #
    #     # Copy files to new directories:
    #     for alignment_directory in alignment_directories:
    #         for alignment_file in glob.glob(f'{alignment_directory}/*.aln.trimmed.fasta'):
    #             shutil.copy(alignment_file, output_alignment_directory)
    #
    #     logger.info(f'{"[INFO]:":10} All trimmed/masked/cut selected alignments have been copied to the directory: '
    #                 f'"{output_alignment_directory}"')
    #
    #     for treefile_directory in treefile_directories:
    #         for treefile in glob.glob(f'{treefile_directory}/*.treefile'):
    #             shutil.copy(treefile, output_treefile_directory)
    #
    #     logger.info(f'{"[INFO]:":10} All trees from trimmed/masked/cut selected alignments have been copied to the '
    #                 f'directory: "{output_treefile_directory}"')
    #
    elif args.from_prune_paralogs:
        collate_pruned_alignments(args.from_prune_paralogs,
                                  directory_contents,
                                  logger=logger)


    # logger.debug(f'alignment_directories are: {alignment_directories}')
    # logger.debug(f'treefile_directories are {treefile_directories}')

    logger.info(f'{"[INFO]:":10} Finished collating alignments and corresponding trees.')
