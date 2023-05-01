# Author: Chris Jackson April 2023 chris.jackson@rbg.vic.gov.au
# https://github.com/chrisjackson-pellicle

"""
Takes a folder of QC'd newick trees and writes the following reports on putative paralogy:

    - Count of duplicated tips with respect to how many tips are present for each locus.
        => A column of loci, a column of # of taxa, and a column of # taxa duplicated.

    - Counts of loci where each taxon is duplicated.
        => A column of taxa, and a column of # of loci where that taxon is duplicated.

"""

import glob
import os
import sys
import textwrap
from collections import defaultdict

from paragone import utils
from paragone import tree_utils
from paragone import newick3


def write_putative_paralogy_reports(treefile_directory,
                                    report_directory,
                                    logger=None):
    """
    Writes per locus and per taxon reports for putative paralogy in the QC'd trees


    :param str treefile_directory: name of input folder containing QC'd tree files
    :param str report_directory: path to directory for report files
    :param logging.Logger logger: a logger object
    :return:
    """

    locus2stats = defaultdict()
    taxon2stats = defaultdict(list)
    all_taxa = set()

    # Read in each tree file:
    for treefile in glob.glob(f'{treefile_directory}/*.treefile'):
        taxon_to_count_dict = defaultdict(list)
        treefile_basename = os.path.basename(treefile)
        tree_name = tree_utils.get_cluster_id(treefile_basename)

        logger.info(f'{"[INFO]:":10} Analysing tree {treefile_basename}...')

        # Read in the tree and check number of taxa:
        with open(treefile, "r") as infile:
            intree = newick3.parse(infile.readline())
            for node in intree.iternodes(order=0):  # PREORDER, root to tip
                if node.istip:
                    base_taxon_name = node.label.split('.')[0]
                    all_taxa.add(base_taxon_name)

                    taxon_to_count_dict[base_taxon_name].append(node.label)

        # Get total number of tips:
        total_tips = [tip_name for tip_name_list in taxon_to_count_dict.values() for tip_name in tip_name_list]
        taxa_with_more_than_one_tip = []
        for taxon, tip_list in taxon_to_count_dict.items():
            if len(tip_list) > 1:
                taxa_with_more_than_one_tip.append(taxon)
                taxon2stats[taxon].append(tree_name)

        locus2stats[tree_name] = len(total_tips), taxa_with_more_than_one_tip

    # Write reports:
    locus_report_outfile = f'{report_directory}/per_locus_paralogy_report_post_tree_qc.tsv'
    taxon_report_outfile = f'{report_directory}/per_taxon_paralogy_report_post_tree_qc.tsv'

    # Per-locus report:
    with open(f'{locus_report_outfile}', 'w') as locus_report_handle:
        locus_report_handle.write(f'locus\tnum_taxa_total\tnum_taxa_>1_tip\t>1_tip_taxa_names\n')
        for locus, stats_tuple in sorted(locus2stats.items()):
            if len(stats_tuple[1]) != 0:
                duplicated_taxa = '; '.join(sorted(stats_tuple[1]))
            else:
                duplicated_taxa = 'None'

            locus_report_handle.write(f'{locus}\t{stats_tuple[0]}\t{len(stats_tuple[1])}\t{duplicated_taxa}\n')

    # Per taxon report:
    with open(f'{taxon_report_outfile}', 'w') as taxon_report_handle:
        taxon_report_handle.write(f'taxon\tloci_where_>1_tip\n')
        for taxon in sorted(list(all_taxa)):
            loci_list = taxon2stats[taxon]
            if len(loci_list) != 0:
                loci_list = '; '.join(sorted(loci_list))
            else:
                loci_list = 'None'
            taxon_report_handle.write(f'{taxon}\t{loci_list}\n')

    return locus_report_outfile, taxon_report_outfile


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

    logger.debug(f'{"[INFO]:":10} Module paralogy_report was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    logger.info(f'\n{"[INFO]:":10} ======> WRITING REPORTS ON PUTATIVE PARALOGY IN QC TREES <======\n')

    treefile_directory = '13_pre_paralog_resolution_trees'
    tree_file_suffix = '.treefile'

    # Checking input directories and files:
    directory_suffix_dict = {treefile_directory: tree_file_suffix}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Write putative paralogy reports:
    locus_report_outfile, taxon_report_outfile = \
        write_putative_paralogy_reports(treefile_directory,
                                        report_directory,
                                        logger=logger)

    fill = textwrap.fill(f'{"[INFO]:":10} Finished writing putative paralogy reports to:'
                         f'\n{locus_report_outfile}\n{taxon_report_outfile}',
                         width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
    logger.info(fill)
