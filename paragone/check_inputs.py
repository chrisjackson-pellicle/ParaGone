#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
- Checks gene names in paralog files and the external outgroup file (if provided) for period/dots, and converts them to
  underscores.
- Checks if there is an outgroup (internal or external) for each gene paralog fasta file, and writes a *.tsv report
  file.
- Writes a  *.tsv file listing all specified 'internal' and 'external' outgroups.
"""

import logging
import sys
import os
import re
import glob
import shutil
from collections import defaultdict
from Bio import SeqIO
import textwrap

from paragone import utils


def sanitise_gene_names(paralogs_folder,
                        gene_name_delimiter,
                        gene_name_field_num,
                        file_of_external_outgroups,
                        logger=None):
    """
    Checks gene names in paralog files and the external outgroup file (if the latter is provided) for periods/dots,
    and converts them to underscores.

    :param str paralogs_folder: path to folder containing input fasta files
    :param str gene_name_delimiter: delimiter in paralog filename to extract gene name. Default is "_".
    :param int gene_name_field_num: from paralog filename, number of fields to extract gene name. Default is 1
    :param str/None file_of_external_outgroups: path to external outgroups file, or None
    :param logging.Logger logger: a logger object
    :return str, str/None sanitised_paralog_output_folder, sanitised_external_outgroup_file,None:
    """

    # Create folder for input fasta files with sanitised names:
    sanitised_input_folder = f'01_input_paralog_fasta_with_sanitised_filenames'
    utils.createfolder(sanitised_input_folder)

    input_fasta_count = 0

    # Sanitise filenames in input paralogs folder:
    fill_1 = textwrap.fill(f'{"[INFO]:":10} Sanitising input paralog fasta filenames. Any dots/periods (".") in the '
                           f'gene name component of the filename will be replaced with underscores ("_"). Sanitised '
                           f'files will be written to directory: "{sanitised_input_folder}".',
                           width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

    fill_2 = textwrap.fill(f'{"[INFO]:":10} Delimiter for gene name extraction: "{gene_name_delimiter}"',
                           width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

    fill_3 = textwrap.fill(f'{"[INFO]:":10} Fields for gene name extraction: {gene_name_field_num}',
                           width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

    logger.info(f'{fill_1}\n{fill_2}\n{fill_3}')

    for file in glob.glob(f'{paralogs_folder}/*.fasta'):
        input_fasta_count += 1
        gene_name_base, ext = os.path.splitext(os.path.basename(file))  # Takes the entire basename as the gene name
        gene_name = '_'.join(gene_name_base.split(gene_name_delimiter)[0:gene_name_field_num])
        gene_name_sanitised = re.sub('[.]', '_', gene_name)
        paralog_filename_sanitised = f'{gene_name_sanitised}{ext}'
        shutil.copy(file, f'{sanitised_input_folder}/{paralog_filename_sanitised}')

    logger.info(f'{"[INFO]:":10} Number of input fasta files: {input_fasta_count}')

    # Sanitise gene names in the external outgroup fasta file, if provided:
    if file_of_external_outgroups:
        sanitised_external_outgroups_filename = f'external_outgroups_sanitised.fasta'

        fill = textwrap.fill(f'{"[INFO]:":10} Sanitising outgroup fasta file. Any dots/periods (".") in the gene name '
                             f'component of the sequence fasta headers will be replaced with underscores ("_").  '
                             f'Sanitised file will be written to: "{sanitised_external_outgroups_filename}"',
                             width=90, subsequent_indent=' ' * 11,  break_on_hyphens=False)

        logger.info(f'{fill}')

        sanitised_seqs_to_write = []
        seqs = SeqIO.parse(file_of_external_outgroups, 'fasta')

        for seq in seqs:
            gene_name = seq.id.split('-')[-1]
            sample_name = seq.id.split('-')[0]
            gene_name_sanitised = re.sub('[.]', '_', gene_name)
            seq_name_sanitised = f'{sample_name}-{gene_name_sanitised}'

            # Re-name sequence:
            seq.id = seq_name_sanitised
            seq.name = seq_name_sanitised
            seq.description = seq_name_sanitised

            sanitised_seqs_to_write.append(seq)

        with open(sanitised_external_outgroups_filename, 'w') as sanitised_outgroups_files:
            SeqIO.write(sanitised_seqs_to_write, sanitised_outgroups_files, 'fasta')

        return sanitised_input_folder, sanitised_external_outgroups_filename

    return sanitised_input_folder, None


def check_outgroup_coverage(folder_of_paralog_files,
                            list_of_internal_outgroups,
                            file_of_external_outgroups,
                            list_of_external_outgroups=None,
                            logger=None):
    """
    Check the number of input genes that have an outgroup sequence in either the list_of_internal_outgroups (i.e.
    corresponding to samples within the existing paralog fasta file), or within a file of external outgroup sequences
    (i.e. new taxa to add as outgroups). Writes a report of gene coverage and the corresponding outgroup(s).

    Writes a log file listing the internal and external outgroup taxa, for using with command `paragone
    align_selected_and_tree`.

    :param str folder_of_paralog_files:
    :param list/None list_of_internal_outgroups:
    :param str/None file_of_external_outgroups:
    :param list/None list_of_external_outgroups:
    :param logging.Logger logger: a logger object
    :return:
    """

    logger.info(f'{"[INFO]:":10} Checking outgroup coverage...')

    logger.debug(f'list_of_internal_outgroups is: {list_of_internal_outgroups}')
    logger.debug(f'list_of_external_outgroups to select from external outgroup file is: {list_of_external_outgroups}')

    # Read in paralog fasta files, and create a dictionary of gene_id:list_of_seq_names:
    paralog_dict = defaultdict(list)
    for fasta in glob.glob(f'{folder_of_paralog_files}/*.fasta'):
        gene_id = os.path.basename(fasta).split('.')[0]  # get prefix e.g. '4471'
        seqs = SeqIO.parse(fasta, 'fasta')
        for seq in seqs:
            paralog_dict[gene_id].append(seq.name)

    # Read in external outgroups file, and create a dictionary of gene_id:list_of_seq_names:
    if file_of_external_outgroups:  # dict not created if no outgroups file provided
        external_outgroup_dict = defaultdict(list)
        seqs = SeqIO.parse(file_of_external_outgroups, 'fasta')
        for seq in seqs:
            gene_id = seq.name.split('-')[-1]  # get gene id e.g. '4471'
            external_outgroup_dict[gene_id].append(seq.name)

    # For each taxon name in list_of_internal_outgroups, check whether there is a sequence for each gene, and create
    # a dictionary:
    internal_outgroup_coverage_dict = defaultdict(list)
    if list_of_internal_outgroups:
        for gene, sequence_names in paralog_dict.items():
            internal_outgroup = 0
            for taxon in list_of_internal_outgroups:
                if taxon in sequence_names or f'{taxon}.main' in sequence_names:  # i.e. if there are paralogs
                    internal_outgroup += 1
                    internal_outgroup_coverage_dict[gene].append(taxon)
                else:
                    logger.debug(f'No internal outgroup sequence for gene {gene} for taxon {taxon}')
            if internal_outgroup == 0:
                internal_outgroup_coverage_dict[gene].append('No internal outgroup')

    # If filtering the external outgroups for specified taxa (i.e. if a list of taxon names to select from the
    # external outgroup file is provided), check whether there is a sequence for each taxon name for each gene,
    # and create a dictionary:
    external_outgroup_coverage_dict = defaultdict(list)
    if file_of_external_outgroups and list_of_external_outgroups:
        for gene, sequences in paralog_dict.items():
            gene_lookup = external_outgroup_dict[gene]  # recovers a list
            if len(gene_lookup) == 0:
                external_outgroup_coverage_dict[gene].append('No external outgroup')
            for taxon in list_of_external_outgroups:
                if taxon in ['-'.join(name.split('-')[:-1]) for name in gene_lookup]:  # i.e the names without the
                    # gene id suffix)
                    external_outgroup_coverage_dict[gene].append(taxon)
                else:
                    logger.debug(f'No external outgroup sequence for gene {gene} for taxon {taxon}')

    if file_of_external_outgroups and not list_of_external_outgroups:  # i.e. if NOT filtering external outgroups
        for gene, sequences in paralog_dict.items():
            gene_lookup = external_outgroup_dict[gene]
            if len(gene_lookup) == 0:
                external_outgroup_coverage_dict[gene].append('No external outgroup')
            for seq in gene_lookup:
                external_outgroup_coverage_dict[gene].append('-'.join(seq.split('-')[:-1]))

    # Iterate over all genes from paralogs dict, and check for internal and/or external outgroups. Write a tsv file
    # of results:
    outgroup_coverage_report = f'00_logs_and_reports/reports/outgroup_coverage_report.tsv'
    outgroup_taxon_list = f'00_logs_and_reports/reports/outgroup_taxon_list.tsv'

    with open(outgroup_coverage_report, 'w') as tsv_report:
        number_of_paralog_files = len(paralog_dict)
        num_paralog_files_with_internal_outgroup = 0
        num_paralog_files_with_external_outgroup = 0
        num_paralog_files_with_both_internal_and_external_outgroup = 0
        tsv_report.write(f'Gene_name\tInternal_outgroup_taxa\tExternal_outgroup_taxa\n')
        for gene, sequences in sorted(paralog_dict.items()):
            internal_outgroups = internal_outgroup_coverage_dict[gene]
            external_outgroups = external_outgroup_coverage_dict[gene]
            if len(internal_outgroups) != 0 and 'No internal outgroup' not in internal_outgroups:
                num_paralog_files_with_internal_outgroup += 1
            if len(external_outgroups) != 0 and 'No external outgroup' not in external_outgroups:
                num_paralog_files_with_external_outgroup += 1
            if len(internal_outgroups) != 0 and 'No internal outgroup' not in internal_outgroups and len(
                    external_outgroups) != 0 and 'No external outgroup' not in external_outgroups:
                num_paralog_files_with_both_internal_and_external_outgroup += 1
            tsv_report.write(f'{gene}\t{";".join(internal_outgroups)}\t{";".join(external_outgroups)}\n')

        fill_1 = textwrap.fill(f'{"[INFO]:":10} Input paralog fasta files with >= one internal outgroup sequence:'
                               f' {num_paralog_files_with_internal_outgroup} of {number_of_paralog_files}',
                               width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        fill_2 = textwrap.fill(f'{"[INFO]:":10} Input paralog fasta files with >= one external outgroup sequence:'
                               f' {num_paralog_files_with_external_outgroup} of {number_of_paralog_files}',
                               width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        fill_3 = textwrap.fill(f'{"[INFO]:":10} Input paralog fasta files with >= one internal AND '
                               f'external outgroup sequence:'
                               f' {num_paralog_files_with_both_internal_and_external_outgroup} of'
                               f' {number_of_paralog_files}',
                               width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
        fill_4 = textwrap.fill(f'{"[INFO]:":10} An outgroup coverage report has been written to file '
                               f'"{outgroup_coverage_report}".',
                               width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

        # Log to stderr, to be captured by Nextflow process to e.g. print a warning
        logger.info(fill_1)
        logger.info(fill_2)
        logger.info(fill_3)
        logger.info(fill_4)

        # Write a log file of 'internal' and 'external' outgroup taxa; used in the command `paragone
        # align_selected_and_tree`:
        internal_ingroups = set(taxon for gene_name, taxon_list in internal_outgroup_coverage_dict.items() for taxon
                                in taxon_list if taxon not in ['No internal outgroup'])

        external_ingroups = set(taxon for gene_name, taxon_list in external_outgroup_coverage_dict.items() for taxon
                                in taxon_list if taxon not in ['No external outgroup'])

        with open(outgroup_taxon_list, 'w') as outgroup_taxon_list_handle:
            for taxon in internal_ingroups:
                outgroup_taxon_list_handle.write(f'INTERNAL_OUTGROUP\t{taxon}\n')
            for taxon in external_ingroups:
                outgroup_taxon_list_handle.write(f'EXTERNAL_OUTGROUP\t{taxon}\n')

        fill = textwrap.fill(f'{"[INFO]:":10} An outgroup taxon list has been written to file "'
                             f'{outgroup_taxon_list}".',
                             width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)

        logger.info(fill)


########################################################################################################################
########################################################################################################################

def main(args, logger=None):
    """
    Entry point for the paragone_main.py script.

    :param args: argparse namespace with subparser options for function main()
    :param logging.Logger logger: a logger object
    :return:
    """

    logger.debug(f'{"[INFO]:":10} Module check_inputs was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)

    logger.debug(f'{fill}')
    logger.debug(args)

    logger.info(f'{"[INFO]:":10} ======> CHECKING INPUT FILES <======\n')

    # Checking input directories and files:
    directory_suffix_dict = {args.gene_fasta_directory: '.fasta'}
    file_list = []

    if args.external_outgroups_file:
        file_list.append(args.external_outgroups_file)

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Check gene names in input paralog files, and the external outgroup file (if provided), for periods/dots,
    # and convert them to underscores:
    paralogs_folder_sanitised, external_outgroups_file_sanitised = \
        sanitise_gene_names(args.gene_fasta_directory,
                            args.gene_name_delimiter,
                            args.gene_name_field_num,
                            args.external_outgroups_file,
                            logger=logger)

    # Check outgroup coverage for each input file:
    check_outgroup_coverage(paralogs_folder_sanitised,
                            args.internal_outgroups,
                            external_outgroups_file_sanitised,
                            list_of_external_outgroups=args.external_outgroups,
                            logger=logger)

    logger.info(f'{"[INFO]:":10} Finished checking input files.\n')
