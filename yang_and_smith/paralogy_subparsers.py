#!/usr/bin/env python

"""
Contains argument subparsers
"""

import textwrap
import logging
import sys


def add_check_and_batch_parser(subparsers):
    """
    Parser for check_and_batch

    :param argparse._SubParsersAction subparsers:
    :return None: no return value specified; default is None
    """

    parser_check_and_batch = subparsers.add_parser('check_and_batch',
                                                   help='Check input file, outgroup coverage, create batches')
    parser_check_and_batch.add_argument('gene_fasta_directory',
                                        type=str,
                                        help="Directory contains fasta files including paralogs")
    parser_check_and_batch.add_argument('--external_outgroups_file',
                                        type=str,
                                        default=None,
                                        help='File in fasta format with additional outgroup sequences to add to each '
                                             'gene')
    parser_check_and_batch.add_argument('--external_outgroup',
                                        action='append',
                                        type=str,
                                        dest='external_outgroups',
                                        default=None,
                                        help='If one or more taxon names are provided, only use these sequences from '
                                             'the user-provided external_outgroups_file')
    parser_check_and_batch.add_argument('--internal_outgroup',
                                        action='append',
                                        type=str,
                                        dest='internal_outgroups',
                                        default=None,
                                        help='Taxon name to use as an internal outgroup (i.e. present in input '
                                             'paralog files')
    parser_check_and_batch.add_argument('--batch_size',
                                        type=int,
                                        default=20,
                                        help='Number of fasta files for each batch of input paralog fasta files. '
                                             'Default is: %(default)s')

    return parser_check_and_batch
