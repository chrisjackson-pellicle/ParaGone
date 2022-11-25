#!/usr/bin/env python

# Author: # Author: Yang and Smith
# Modified by: Chris Jackson chris.jackson@rbg.vic.gov.au

"""
If no outgroup, only output ingroup clades with no taxon repeats
If outgroup present, extract rooted ingroup clades and prune paralogs

Prepare a taxon file, with each line look like (separated by tab):
IN	taxonID1
IN	taxonID2
OUT	taxonID3
"""

import glob
import os
import sys
import textwrap
from collections import defaultdict

from paragone import tree_utils
from paragone import newick3
from paragone import utils


def prune_paralogs_rt(curroot,
                      treefile_basename,
                      output_file_id,
                      ingroups,
                      outgroups,
                      min_ingroup_taxa,
                      logger=None):
    """

    :param curroot:
    :param treefile_basename:
    :param output_file_id:
    :param ingroups:
    :param outgroups:
    :param min_ingroup_taxa:
    :param logging.Logger logger: a logger object
    :return collections.defaultdict inclades_with_fewer_than_min_ingroup_taxa: a dictionary of treename:[inclades
    with fewer than min_ingroup_taxa]
    """

    # Recover a list of rooted inclades:
    inclades_list, inclades_with_fewer_than_min_ingroup_taxa = \
        tree_utils.extract_rooted_ingroup_clades(curroot,
                                                 treefile_basename,
                                                 ingroups,
                                                 outgroups,
                                                 min_ingroup_taxa,
                                                 logger=logger)

    ortho_dict = defaultdict(list)
    ortho_dict_fewer_than_min_taxa = defaultdict(list)

    inclade_count = 0

    for inclade in inclades_list:
        inclade_count += 1

        inclade_name = f'{output_file_id}.inclade{str(inclade_count)}'
        with open(inclade_name, "w") as outfile:
            outfile.write(newick3.tostring(inclade) + ";\n")

        # Do something here:
        orthologs = tree_utils.get_ortho_from_rooted_inclade(inclade,
                                                             logger=logger)

        ortho_count = 0
        for ortho in orthologs:
            if len(tree_utils.get_front_labels(ortho)) >= min_ingroup_taxa:
                ortho_count += 1

                ortho_name = f'{inclade_name}.ortho{str(ortho_count)}.tre'
                with open(ortho_name, "w") as outfile:
                    outfile.write(newick3.tostring(ortho) + ";\n")
                ortho_dict[f'inclade_{inclade_count}'].append(newick3.tostring(ortho))
            else:
                logger.debug(f'Ortho from rooted ingroup clade for tree {treefile_basename} has fewer taxa than the '
                             f'min_ingroup_taxa os {min_ingroup_taxa}. Skipping ortho...')

                ortho_dict_fewer_than_min_taxa[f'inclade_{inclade_count}'].append(newick3.tostring(ortho))

    return inclades_list, inclades_with_fewer_than_min_ingroup_taxa, ortho_dict, ortho_dict_fewer_than_min_taxa


def write_rt_report(treefile_directory,
                    tree_stats_collated_dict,
                    logger=None):
    """
    Writes a *.tsv report detailing for RooTed outgroup (RT) pruning process.

    :param str treefile_directory: name of tree file directory for report filename
    :param tree_stats_collated_dict:
    :param logging.Logger logger: a logger object
    :return:
    """

    basename = os.path.basename(treefile_directory)
    report_filename = f'00_logs_and_reports_resolve_paralogs/reports/{basename.lstrip("17_")}_RT_report.tsv'

    logger.info(f'{"[INFO]:":10} Writing RooTed outgroup (RT) report to file {report_filename}')

    trees_with_unrecognised_names_count = 0
    trees_with_fewer_than_min_ingroup_taxa_count = 0
    trees_with_no_outgroup_but_1to1_orthologs_count = 0
    trees_with_no_outgroup_and_putative_paralogs_count = 0
    trees_with_inclades_above_min_taxa_taxa_count = 0
    trees_with_inclades_below_min_taxa_taxa_count = 0
    trees_with_inclade_to_ortho_above_min_taxa_count = 0
    trees_with_inclade_to_ortho_below_min_taxa_count = 0

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
            check = dictionaries['no_outgroup_taxa_1to1_orthologs']
            trees_with_no_outgroup_but_1to1_orthologs_count += 1
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries['no_outgroup_taxa_putative_paralogs']
            trees_with_no_outgroup_and_putative_paralogs_count += 1
            tree_stats.append('Y')
        except KeyError:
            tree_stats.append('N')

        try:
            check = dictionaries['inclades_above_min_taxa']
            assert len(check) > 0
            trees_with_inclades_above_min_taxa_taxa_count += 1
            tree_stats.append(len(check))
        except AssertionError:
            tree_stats.append('0')
        except KeyError:
            tree_stats.append('N/A')

        try:
            check = dictionaries['inclades_below_min_taxa']
            assert len(check) > 0
            trees_with_inclades_below_min_taxa_taxa_count += 1
            tree_stats.append(len(check))
        except AssertionError:
            tree_stats.append('0')
        except KeyError:
            tree_stats.append('N/A')

        try:
            check = dictionaries['inclade_to_ortho_above_min_taxa_dict']
            assert len(check) > 0
            ortho_count = 0
            for key, values, in check.items():
                ortho_count += len(values)
            trees_with_inclade_to_ortho_above_min_taxa_count += 1
            tree_stats.append(ortho_count)
        except AssertionError:
            tree_stats.append('0')
        except KeyError:
            tree_stats.append('N/A')

        try:
            check = dictionaries['inclade_to_ortho_below_min_taxa_dict']
            assert len(check) > 0
            ortho_count = 0
            for key, values, in check.items():
                ortho_count += len(values)
            trees_with_inclade_to_ortho_below_min_taxa_count += 1
            tree_stats.append(ortho_count)
        except AssertionError:
            tree_stats.append('0')
        except KeyError:
            tree_stats.append('N/A')

        all_tree_stats_for_report.append(tree_stats)

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'\t'
                            f'Unrecognised taxa tree (skipped)\t'
                            f'< than minimum ingroup taxa (tree skipped)\t'
                            f'No outgroup taxa and putative paralogs (tree skipped)\t'
                            f'No outgroup taxa but 1to1 orthologs\t'
                            f'Rooted inclades recovered\t'
                            f'Rooted inclades < minimum ingroup taxa\t'
                            f'Ortholog groups recovered from rooted inclades > minimum ingroup taxa\t'
                            f'Ortholog groups recovered from rooted inclades < minimum taxa'
                            f'\n')

        report_handle.write(f'Number of trees\t'
                            f'{trees_with_unrecognised_names_count}\t'
                            f'{trees_with_fewer_than_min_ingroup_taxa_count}\t'
                            f'{trees_with_no_outgroup_but_1to1_orthologs_count}\t'
                            f'{trees_with_no_outgroup_and_putative_paralogs_count}\t'
                            f'{trees_with_inclades_above_min_taxa_taxa_count}\t'
                            f'{trees_with_inclades_below_min_taxa_taxa_count}\t'
                            f'{trees_with_inclade_to_ortho_above_min_taxa_count}\t'
                            f'{trees_with_inclade_to_ortho_below_min_taxa_count}'
                            f'\n')

        for stats in all_tree_stats_for_report:
            stats_joined = '\t'.join([str(stat) for stat in stats])
            report_handle.write(f'{stats_joined}\n')


def main(args):
    """
    Entry point for the paragone.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = utils.setup_logger(__name__, '00_logs_and_reports_resolve_paralogs/logs/12_prune_paralogs_RT')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand prune_paralogs_rt was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')
    logger.debug(args)

    # Checking input directories and files:
    directory_suffix_dict = {args.treefile_directory: args.tree_file_suffix}
    file_list = []

    if args.in_and_outgroup_list:
        file_list.append(args.in_and_outgroup_list)

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Create output folder for pruned trees:
    output_folder = f'19_{os.path.basename(args.treefile_directory).lstrip("17_")}_pruned_RT'
    utils.createfolder(output_folder)

    # Parse the ingroup and outgroup text file:
    ingroups, outgroups = utils.parse_ingroup_and_outgroup_file(args.in_and_outgroup_list,
                                                                logger=logger)

    # Create dict for report file:
    tree_stats_collated = defaultdict(lambda: defaultdict())

    # Iterate over tree and prune with RT algorithm:
    for treefile in glob.glob(f'{args.treefile_directory}/*{args.tree_file_suffix}'):
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
                logger.warning(f'{"[WARNING]:":10} Taxon names {unrecognised_names} in tree {treefile_basename} not '
                               f'found in ingroups or outgroups. Skipping tree...')

                tree_stats_collated[treefile_basename]['unrecognised_names'] = unrecognised_names
                continue

            # Check if tree contains more than the minimum number of taxa, else skip:
            if len(ingroup_names) < args.minimum_taxa:
                logger.warning(f'{"[WARNING]:":10} Tree {treefile_basename} contains {len(ingroup_names)} taxa; '
                               f'minimum_taxa required is {args.minimum_taxa}. Skipping tree...')

                tree_stats_collated[treefile_basename]['fewer_than_min_ingroup_taxa'] = newick3.tostring(curroot)
                continue

            # If no outgroup at all, but no putative paralogs, write unrooted tree to file:
            if len(outgroup_names) == 0 and len(names) == num_taxa:
                outfile_filename = f'{output_file_id}.unrooted-ortho.tre'

                logger.info(f'{"[INFO]:":10} Tree {treefile_basename} contains no outgroup taxa, but no putative '
                            f'paralogs detected. Writing tree to {outfile_filename}')
                with open(f'{outfile_filename}', 'w') as outfile:
                    outfile.write(newick3.tostring(curroot) + ";\n")

                tree_stats_collated[treefile_basename]['no_outgroup_taxa_1to1_orthologs'] = newick3.tostring(curroot)

            # If no outgroup at all, and putative paralogs detected, skip tree:
            elif len(outgroup_names) == 0:
                logger.info(f'{"[WARNING]:":10} Tree {treefile_basename} contains no outgroup taxa and putative '
                            f'paralogs detected. Skipping tree...')

                tree_stats_collated[treefile_basename]['no_outgroup_taxa_putative_paralogs'] = newick3.tostring(curroot)
                continue

            # If outgroups, run RooTed outgroups (RT) algorithm:
            else:
                inclades_list, \
                inclades_with_fewer_than_min_ingroup_taxa_list, \
                ortho_dict, \
                ortho_dict_fewer_than_min_taxa = \
                    prune_paralogs_rt(curroot,
                                      treefile_basename,
                                      output_file_id,
                                      ingroups,
                                      outgroups,
                                      args.minimum_taxa,
                                      logger=logger)

                # Collate per tree data in dictionary:
                tree_stats_collated[treefile_basename]['inclades_above_min_taxa'] = inclades_list
                tree_stats_collated[treefile_basename]['inclades_below_min_taxa'] = \
                    inclades_with_fewer_than_min_ingroup_taxa_list
                tree_stats_collated[treefile_basename]['inclade_to_ortho_above_min_taxa_dict'] = ortho_dict
                tree_stats_collated[treefile_basename]['inclade_to_ortho_below_min_taxa_dict'] = \
                    ortho_dict_fewer_than_min_taxa

    write_rt_report(args.treefile_directory,
                    tree_stats_collated,
                    logger=logger)

    logger.info(f'{"[INFO]:":10} Finished extracting putative ortholog trees using the RooTed outgroups (RT) '
                f'algorithm.')
