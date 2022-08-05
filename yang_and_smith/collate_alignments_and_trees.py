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

from yang_and_smith import utils


def main(args):
    """
    Entry point for the resolve_paralogs.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = utils.setup_logger(__name__, '00_logs_and_reports_resolve_paralogs/logs/04_collate_alignments_and_trees')

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

    # Create output folders:
    output_alignment_directory = f'06_all_hmmcleaned_alignments'
    utils.createfolder(output_alignment_directory)

    output_treefile_directory = f'07_all_hmmcleaned_alignment_trees'
    utils.createfolder(output_treefile_directory)

    directory_contents = os.listdir('.')
    alignment_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
                             directory.startswith('04_batch')]
    treefile_directories = [directory for directory in directory_contents if os.path.isdir(directory) and
                            directory.startswith('05_batch')]

    logger.debug(f'alignment_directories are: {alignment_directories}')
    logger.debug(f'treefile_directories are {treefile_directories}')

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

    logger.info(f'{"[INFO]:":10} Finished collating cleaned alignments and corresponding trees.')
