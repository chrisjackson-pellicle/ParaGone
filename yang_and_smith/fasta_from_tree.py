# Author: Alexander Schmidt-Lebuhn (modified by Chris Jackson July 2021 chris.jackson@rbg.vic.gov.au)

"""
Takes a newick tree and a fasta aligmment as input. Returns a fasta alignment containing only the sequences
corresponding to the tree tip names.
"""

from Bio import AlignIO
from Bio import Phylo
import Bio.Align
import glob
import os
import sys
import textwrap
import shutil

from yang_and_smith import utils


def subsample_alignments(treefile_directory,
                         tree_suffix,
                         alignment_directory,
                         from_cut_deep_paralogs=False,
                         logger=None):
    """
    Takes a pruned/QC'd tree file, finds the original matching alignment, and sub-samples that alignment to recover
    only sequences corresponding to tree tip names.

    :param str treefile_directory: path to directory containing tree newick files
    :param str tree_suffix: suffix for the tree files
    :param str alignment_directory: path to the directory containing fasta alignments
    :param bool from_cut_deep_paralogs: if True, process tree file names accordingly to recover gene names
    :param logging.Logger logger: a logger object
    :return str, dict output_folder, alignment_filtering_dict: path the output folder with filtered alignments,
    dictionary of filtering stats for each tree/alignment
    """

    logger.info(f'{"[INFO]:":10} Recovering alignment sequences corresponding to tree tip names...')

    treefile_directory_basename = os.path.basename(treefile_directory)
    output_folder = f'{treefile_directory_basename}_alignments'
    utils.createfolder(output_folder)

    # Capture number of sequences pre and post filtering in a dictionary for report:
    alignment_filtering_dict = {}

    for tree in glob.glob(f'{treefile_directory}/*{tree_suffix}'):
        read_tree = Phylo.read(tree, "newick")
        tree_terminals = read_tree.get_terminals()
        # print(len(tree_terminals))
        tree_basename = os.path.basename(tree)

        # Derive the matching alignment file name depending on input tree file name:
        if from_cut_deep_paralogs:  # e.g. 4471_1.subtree
            alignment_prefix = '_'.join(tree_basename.split('_')[0:-1])
            output_alignment_prefix = tree_basename.split('.')[0]
            # print(f'alignment_prefix is: {alignment_prefix}')
            matching_alignment = f'{alignment_directory}/{alignment_prefix}.paralogs.aln.hmm.trimmed.fasta'
            # print(matching_alignment)
        else:  # e.g. 4691_1.1to1ortho.tre, 4471_1.inclade1.ortho1.tre, 4527_1.MIortho1.tre. etc
            alignment_prefix = tree_basename.split('.')[0]
            output_alignment_prefix = '.'.join(tree_basename.split('.')[0:-1])
            # print(f'alignment_prefix is: {alignment_prefix}')
            matching_alignment = f'{alignment_directory}/{alignment_prefix}.outgroup_added.aln.trimmed.fasta'
            # print(f'matching_alignment is: {matching_alignment}')

        # Read in original alignments and select seqs matching tree termini:
        alignment = AlignIO.read(matching_alignment, "fasta")
        subalignment = Bio.Align.MultipleSeqAlignment([])
        for k in range(0, len(alignment)):
            for j in range(0, len(tree_terminals)):
                if tree_terminals[j].name == alignment[k].id:
                    subalignment.append(alignment[k])

        # Capture data:
        alignment_filtering_dict[tree_basename] = [len(tree_terminals), len(alignment), len(subalignment)]

        # Write an alignment of the sub-selected sequences:
        AlignIO.write(subalignment, f'{output_folder}/{output_alignment_prefix}.selected.fasta', "fasta")

    return output_folder, alignment_filtering_dict


def batch_input_files(gene_fasta_directory,
                      output_directory,
                      batch_size=20,
                      logger=None):
    """
    Takes a folder of fasta files, and splits them in to batch folders according to the number provided by
    parameter batch_size.

    :param str gene_fasta_directory: path to input fasta files with sanitised filenames
    :param str output_directory: name of output directory to create
    :param int batch_size: number of fasta files per batch; default is 20
    :param logging.Logger logger: a logger object
    :return:
    """

    utils.createfolder(output_directory)

    fasta_file_list = glob.glob(f'{gene_fasta_directory}/*.fasta')

    def chunks(lst, n):
        """Yield successive n-sized chunks from lst."""
        for i in range(0, len(lst), n):
            yield lst[i:i + n]

    batches = list(chunks(fasta_file_list, batch_size))
    # logger.info(f'Batches are: {batches}')
    batch_num = 1
    for batch in batches:
        utils.createfolder(f'{output_directory}/batch_{batch_num}')
        for fasta_file in batch:
            shutil.copy(fasta_file, f'{output_directory}/batch_{batch_num}')
        batch_num += 1


def write_fasta_from_tree_report(alignment_filtering_dict,
                                 treefile_directory,
                                 from_cut_deep_paralogs,
                                 logger=None):
    """
    Writes a *.tsv report detailing number of tips in QC'd tree, number of sequences in original QC'd alignment,
    and number of sequences in filtered alignment.

    :param dict alignment_filtering_dict: dictionary of filtering stats for each tree/alignment
    :param str treefile_directory: name of tree file directory for report filename
    :param bool from_cut_deep_paralogs: if True, add 'cut' to report filename
    :param logging.Logger logger: a logger object
    :return:
    """

    basename = os.path.basename(treefile_directory)
    if from_cut_deep_paralogs:
        report_filename = f'{basename}_fasta_from_tree_report_cut.tsv'
    else:
        report_filename = f'{basename}_fasta_from_tree_report_final.tsv'

    logger.info(f'{"[INFO]:":10} Writing fasta from tree report to file {report_filename}')

    with open(report_filename, 'w') as report_handle:
        report_handle.write(f'QC tree file\tNumber of tree tips\tNumber seqs in original alignment\tNumber seqs '
                            f'filtered alignment\n')

        for tree_name, stats in alignment_filtering_dict.items():
            report_handle.write(f'{tree_name}\t{stats[0]}\t{stats[1]}\t{stats[2]}\n')


def main(args):
    """
    Entry point for the resolve_paralogs.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    if args.from_cut_deep_paralogs:
        logger = utils.setup_logger(__name__, 'logs_resolve_paralogs/07_fasta_from_tree_cut')
    else:
        logger = utils.setup_logger(__name__, 'logs_resolve_paralogs/07_fasta_from_tree_final')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand fasta_from_tree was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')
    logger.debug(args)

    # Create output folder:
    treefile_directory_basename = os.path.basename(args.treefile_directory)
    output_folder = f'{treefile_directory_basename}_alignments_batches'
    utils.createfolder(output_folder)

    # Recover alignments with sequences corresponding to tree tip names:
    filtered_alignments_folder,  alignment_filtering_dict = \
        subsample_alignments(args.treefile_directory,
                             args.tree_file_suffix,
                             args.alignment_directory,
                             from_cut_deep_paralogs=args.from_cut_deep_paralogs,
                             logger=logger)

    # Batch fasta files for alignment and tree-building steps:
    batch_input_files(filtered_alignments_folder,
                      output_folder,
                      batch_size=args.batch_size,
                      logger=logger)

    # Write a report of pre-and-post filtering stats for each tree/alignments:
    write_fasta_from_tree_report(alignment_filtering_dict,
                                 args.treefile_directory,
                                 args.from_cut_deep_paralogs,
                                 logger=logger)
