#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
Input is fasta files that have undergone QC processes; this QC may have removed internal outgroups (if specified).
This script:

    - Checks for any paralogs in internal outgroups (if the latter are specified), and selects a single
      representative sequence for each taxon  - this will be the sequence with the highest distance from a sample of up
      to 10 ingroup sequences.
    - Aligns the sequences using MAFFT.
    - Generates trees from the alignments using either IQTREE or FastTreeMP.

NOTE: trees generated with FASTREEMP can contain polytomies, which will need to be resolved before processing with
the Yang and Smith scripts!

"""

import sys
import os
import re
import textwrap
import glob
import subprocess
import copy
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Align import MultipleSeqAlignment
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait, as_completed
import traceback

from paragone import utils
from paragone.align_and_clean import run_trimal


def add_outgroup_seqs(qc_alignment_directory,
                      selected_alignment_directory,
                      list_of_internal_outgroups,
                      list_of_external_outgroups,
                      logger=None):
    """
    Check the number of genes that have an outgroup sequence in either the list_of_internal_outgroups (i.e.
    corresponding to samples within the existing [ORIGINAL - might have been pruned out by this stage] paralog fasta
    file), or within a file of external outgroup sequences (i.e. new taxa to add as outgroups). Add internal (if
    removed by QC steps) and external outgroup seqs to each alignment

    Write a tab-separated outgroup file for the Y&S pruning scripts, of the form:

        IN  Euchiton_limosus
        IN  Euchiton_sphaericus
        IN  Pterochaeta_paniculata
        OUT sunf
        etc...

    :param str qc_alignment_directory:
    :param str selected_alignment_directory:
    :param list list_of_internal_outgroups:
    :param list list_of_external_outgroups:
    :param logging.Logger logger: a logger object
    :return str output_folder: name of output folder containing fasta with outgroups added
    """

    logger.debug(f'list_of_internal_outgroups: {list_of_internal_outgroups}')
    logger.debug(f'list_of_external_outgroups_to_select: {list_of_external_outgroups}')

    if len(list_of_internal_outgroups) == 0 and len(list_of_external_outgroups) == 0:
        logger.warning(f'{"[WARNING]:":10} No external or internal outgroups supplied!')

    output_folder = f'10_sequences_from_qc_outgroups_added'
    output_folder = utils.createfolder(output_folder)  # for the outgroups added fasta files

    # Read in original paralog fasta files, and create a dictionary of gene_id:list_of_seq_names for taxa in
    # list_of_internal_outgroups:
    internal_outgroup_dict = defaultdict(lambda: defaultdict(list))
    all_paralog_taxon_names = set()
    internal_outgroup_found_dict = dict()

    for original_alignment in glob.glob(f'{qc_alignment_directory}/*.fasta'):
        gene_id = os.path.basename(original_alignment).split('.')[0]  # get prefix e.g. '4471'
        seqs = SeqIO.parse(original_alignment, 'fasta')
        internal_outgroup_found_dict[gene_id] = False  # set default

        # Populate the internal_outgroup_dict (potentially more than one sequence per taxon):
        for seq in seqs:
            seq_name_prefix = seq.name.split('.')[0]  # assumes there are no other dots ('.') in the sequence name
            all_paralog_taxon_names.add(seq_name_prefix)
            if list_of_internal_outgroups and seq_name_prefix in list_of_internal_outgroups:
                internal_outgroup_dict[gene_id][seq_name_prefix].append(seq)
                internal_outgroup_found_dict[gene_id] = True

        if list_of_internal_outgroups:
            alignment = AlignIO.read(original_alignment, 'fasta')

            # Create a MultipleSeqAlignment object for ingroup sequences only, and take the first 10 sequences from
            # the ingroup alignment for distance matrix calculations:
            alignment_ingroup_seqs_only = []
            ingroup_count = 0
            for sequence in alignment:
                gene_name = sequence.id.split('.')[0]
                if gene_name not in list_of_internal_outgroups:
                    ingroup_count += 1
                    if ingroup_count <= 10:
                        alignment_ingroup_seqs_only.append(sequence)
                    else:
                        break

            alignment_ingroup_only = MultipleSeqAlignment(alignment_ingroup_seqs_only)

            # Filter internal outgroups with more than one seq per taxon (putative paralogs/alleles) to retain the
            # single sequence most divergent from the ingroup sample:
            logger.info(f'{"[INFO]:":10} Detecting paralogs/alleles in internal outgroup sequences for gene'
                        f' {gene_id}...')
            internal_outgroup_dict_filtered = filter_internal_outgroups(internal_outgroup_dict,
                                                                        alignment_ingroup_only,
                                                                        gene_id,
                                                                        logger=logger)
            # Reassign internal_outgroup_dict[gene_id] dictionary entry to selected sequences:
            internal_outgroup_dict[gene_id] = internal_outgroup_dict_filtered[gene_id]

    # Check that internal outgroups were actually found in the "qc_alignment_directory" provided:
    genes_with_internal_outgroups_found = [gene_id for gene_id in internal_outgroup_found_dict.keys() if
                                           internal_outgroup_found_dict[gene_id]]
    genes_with_no_internal_outgroups_found = [gene_id for gene_id in internal_outgroup_found_dict.keys() if not
                                              internal_outgroup_found_dict[gene_id]]

    logger.debug(f'genes_with_internal_outgroups_found: {genes_with_internal_outgroups_found}')
    logger.debug(f'genes_with_no_internal_outgroups_found: {genes_with_no_internal_outgroups_found}')

    if len(list_of_internal_outgroups) != 0 and len(genes_with_internal_outgroups_found) == 0:
        fill = textwrap.fill(f'{"[ERROR]:":10} No internal outgroup sequences were found for any gene within the '
                             f'directory provided: "{qc_alignment_directory}". Please check you have provided the '
                             f'correct directory!',
                             width=90, subsequent_indent=' ' * 11,
                             break_on_hyphens=False)
        logger.error(f'{fill}')
        sys.exit()

    if len(genes_with_no_internal_outgroups_found) != 0:
        fill = textwrap.fill(f'{"[INFO]:":10} No internal outgroup sequences were found for some genes in the '
                             f'directory "{qc_alignment_directory}". This may be expected i.e. if some of your input '
                             f'paralog fasta files do not have any internal outgroup sequences. The genes are:',
                             width=90, subsequent_indent=' ' * 11,
                             break_on_hyphens=False)

        logger.info(f'{fill}')

        for gene_id in genes_with_no_internal_outgroups_found:
            logger.info(f'{" " * 11}{gene_id}')

    # Read in external outgroups file if present, and create a dictionary of gene_id:list_of_seq_names, either for all
    # seqs if no external outgroup taxa specified, or for specified taxa only:
    external_outgroup_dict = defaultdict(list)

    # if file_of_external_outgroups:  # dict not populated if no outgroups file provided
    if list_of_external_outgroups != 0:
        try:
            seqs = SeqIO.parse('external_outgroups_sanitised.fasta', 'fasta')
            for seq in seqs:
                gene_id = seq.name.split('-')[-1]  # get gene id e.g. '4471'
                taxon = '-'.join(seq.name.split('-')[:-1])  # e.g. 'AMBTR'
                # if list_of_external_outgroups:
                if taxon in list_of_external_outgroups:
                    seq.name = taxon  # i.e. we don't want the suffix e.g. '4471' present in the fasta file
                    seq.id = taxon
                    external_outgroup_dict[gene_id].append(seq)
                # else:
                #     seq.name = taxon  # i.e. we don't want the suffix e.g. '4471' present in the fasta file
                #     seq.id = taxon
                #     external_outgroup_dict[gene_id].append(seq)
        except FileNotFoundError:
            logger.error(f'Expected file "external_outgroups_sanitised.fasta" not found!')

    all_external_outgroup_taxon_names = \
        set([seq.name for gene_id, seq_list in external_outgroup_dict.items() for seq in seq_list])

    # Read in QC'd paralog files, add outgroup seqs, and write new fasta files ready for alignment:
    for selected_alignment in glob.glob(f'{selected_alignment_directory}/*.selected.fasta'):
        gene_id = re.sub('_[1-9]$', '', str(os.path.basename(selected_alignment).split('.')[0]))  # 4527_1.selected -> 4527
        gene_id_with_subtree_number = os.path.basename(selected_alignment).split('.')[0]  # 4527_1.selected -> 4527_1

        seqs = list(SeqIO.parse(selected_alignment, 'fasta'))
        if list_of_internal_outgroups:
            # Get seqs that are NOT internal outgroups:
            seqs_to_write = [seq for seq in seqs if seq.name.split('.')[0] not in list_of_internal_outgroups]
        else:
            seqs_to_write = seqs

        for taxon, sequence_list in internal_outgroup_dict[gene_id].items():  # add internal outgroup seqs
            seqs_to_write.extend(sequence_list)

        external_outgroup_seqs = external_outgroup_dict[gene_id]
        seqs_to_write.extend(external_outgroup_seqs)  # add external outgroup seqs

        # Write new files with outgroup sequences added:
        with open(f'{output_folder}/{gene_id_with_subtree_number}.outgroup_added.fasta', 'w') as outgroup_added:
            SeqIO.write(seqs_to_write, outgroup_added, 'fasta')

    # Write the IN and OUT taxon text file required by some paralogy resolution methods (MO, RT):
    if list_of_internal_outgroups:
        ingroup_taxon_names = [name for name in all_paralog_taxon_names if name not in list_of_internal_outgroups]
        logger.debug(f'ingroup_taxon_names: {ingroup_taxon_names}')
    else:
        ingroup_taxon_names = [name for name in all_paralog_taxon_names]
        logger.debug(f'ingroup_taxon_names: {ingroup_taxon_names}')
    with open(f'00_logs_and_reports/reports/in_and_outgroups_list.tsv', 'w') as group_list:
        if list_of_internal_outgroups:
            for taxon in list_of_internal_outgroups:
                group_list.write(f'OUT\t{taxon}\n')
        for taxon in all_external_outgroup_taxon_names:
            group_list.write(f'OUT\t{taxon}\n')
        for taxon in ingroup_taxon_names:
            group_list.write(f'IN\t{taxon}\n')

    return output_folder


def filter_internal_outgroups(internal_outgroup_dict,
                              alignment_ingroup_only,
                              gene_id,
                              logger=None):
    """
    Filters internal outgroups with paralogs to retain the single sequence most divergent from the ingroup

    :param internal_outgroup_dict:
    :param alignment_ingroup_only:
    :param gene_id:
    :param logging.Logger logger: a logger object
    :return:
    """

    # Make a copy of the internal_outgroup_dict so that items can be deleted and the dict returned:
    internal_outgroup_dict_copy = copy.deepcopy(internal_outgroup_dict)

    # Check of there are paralogs for any internal outgroup:
    for taxon_id, sequence_list in internal_outgroup_dict[gene_id].items():
        if len(sequence_list) > 1:
            fill = textwrap.fill(f'{"[INFO]:":10} Taxon {taxon_id} is an internal outgroup, and has more than one '
                                 f'sequence for gene {gene_id}. Only the sequence most divergent from the ingroup '
                                 f'taxa will be retained.',
                                 width=90, subsequent_indent=' ' * 11, break_on_hyphens=False)
            logger.info(f'{fill}')

            seq_to_keep_distance = 0.0  # default value, same at every position
            seq_to_keep = None

            for seq in sequence_list:
                alignment_ingroup_only_edit = copy.deepcopy(alignment_ingroup_only)
                assert len(seq) == alignment_ingroup_only.get_alignment_length()
                alignment_ingroup_only_edit.append(seq)  # add single sequence to ingroup alignment

                for sequence in alignment_ingroup_only_edit:  # Distance matrix requires capital letters
                    sequence.seq = sequence.seq.upper()

                # Create a distance matrix
                skip_letters = set(letter for sequence in alignment_ingroup_only_edit for letter in sequence.seq if
                                   letter not in ['A', 'T', 'C', 'G'])
                my_calculator = DistanceCalculator('blastn', skip_letters=''.join(skip_letters))
                trimmed_dm = my_calculator.get_distance(alignment_ingroup_only_edit)
                distance_values = trimmed_dm[seq.name]
                sorted_distance_values = sorted(distance_values, key=float)
                closest_sequence_distance = sorted_distance_values[1]  # skip zero as it's self-vs-self

                # Make sure a sequence is selected:
                if not seq_to_keep:
                    seq_to_keep_distance = closest_sequence_distance
                    seq_to_keep = seq
                    continue

                # If a sequence has been selected already, compare it to the current seq:
                if closest_sequence_distance > seq_to_keep_distance:  # i.e. it's more diverged
                    logger.debug(f'{closest_sequence_distance} is larger than {seq_to_keep_distance}')
                    seq_to_keep_distance = closest_sequence_distance
                    seq_to_keep = seq

            logger.debug(f'Keeping sequence {seq_to_keep.name} with distance {seq_to_keep_distance}')
            internal_outgroup_dict_copy[gene_id][taxon_id] = [seq_to_keep]

    return internal_outgroup_dict_copy


def mafft_align_multiprocessing(fasta_to_align_folder,
                                algorithm='auto',
                                pool_threads=1,
                                mafft_threads=1,
                                logger=None):

    """
    Generate alignments via function <mafft_align> using multiprocessing.

    :param str fasta_to_align_folder: path to folder containing input fasta files with outgroups added
    :param str algorithm: algorithm to use for mafft alignment; default is 'auto'
    :param int pool_threads: number of alignments to run concurrently
    :param int mafft_threads: number of threads to use for each concurrent alignment
    :param logging.Logger logger: a logger object
    :return str output_folder: name of the output folder containing alignments
    """

    output_folder = f'11_pre_paralog_resolution_alignments'
    utils.createfolder(output_folder)

    logger.info(f'{"[INFO]:":10} Generating alignments from fasta files using MAFFT...')

    # Filter out any input files with fewer than four sequences:
    target_genes = []
    for fasta_file in sorted(glob.glob(f'{fasta_to_align_folder}/*.outgroup_added.fasta')):
        with open(fasta_file, 'r') as input_fasta_handle:
            seqs = list(SeqIO.parse(input_fasta_handle, 'fasta'))
            if len(seqs) < 4:
                logger.warning(f'{"[WARNING]:":10} Skipping file {fasta_file} as it contains fewer than 4 sequences!')
                continue
            else:
                target_genes.append(fasta_file)

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(mafft_align,
                                      fasta_file,
                                      algorithm,
                                      output_folder,
                                      counter,
                                      lock,
                                      num_files_to_process=len(target_genes),
                                      threads=mafft_threads,
                                      logger=logger)
                          for fasta_file in target_genes]

        for future in as_completed(future_results):

            try:
                check = future.result()

            except Exception as error:
                print(f'Error raised: {error}')
                tb = traceback.format_exc()
                print(f'traceback is:\n{tb}')

        wait(future_results, return_when="ALL_COMPLETED")

    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                      utils.file_exists_and_not_empty(alignment)]

    logger.debug(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')

    return output_folder


def mafft_align(fasta_file,
                algorithm,
                output_folder,
                counter,
                lock,
                num_files_to_process,
                threads=1,
                logger=None):
    """
    Use mafft to align a fasta file of sequences, using the algorithm (default is 'auto') and number of threads
    provided. Returns filename of the output alignment.

    :param str fasta_file: path to a fasta file
    :param str algorithm: algorithm to use for mafft alignment; default is 'auto'
    :param str output_folder: name of output folder for alignments
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int num_files_to_process: total number of fasta files for alignment
    :param int threads: number of threads to use for alignment program
    :param logging.Logger logger: a logger object
    :return str expected_alignment_file/expected_alignment_file_trimmed: filename of output alignment
    """

    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'

    try:
        assert utils.file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
        with lock:
            counter.value += 1

        return os.path.basename(expected_alignment_file)

    except AssertionError:

        if algorithm == 'auto':
            command = f'mafft --auto --thread {threads} {fasta_file} > {expected_alignment_file}'
        else:
            command = f'{algorithm} --thread {threads} {fasta_file} > {expected_alignment_file}'

        logger.debug(f'{"[INFO]:":10} Performing MAFFT alignment with command: {command}')

        try:

            result = subprocess.run(command,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    check=True,
                                    shell=True)

            logger.debug(f'MAFFT check_returncode() is: {result.check_returncode()}')
            logger.debug(f'MAFFT stdout is: {result.stdout}')
            logger.debug(f'MAFFT stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'MAFFT FAILED. Output is: {exc}')
            logger.error(f'MAFFT stdout is: {exc.stdout}')
            logger.error(f'MAFFT stderr is: {exc.stderr}')

            raise ValueError('There was an issue running MAFFT. Check input files!')

        with lock:
            counter.value += 1
            logger.debug(f'Aligned file {fasta_file_basename}')

        return os.path.basename(expected_alignment_file)

    finally:
        with lock:
            sys.stderr.write(f'\r{"[INFO]:":10} Finished generating alignment for file {fasta_file_basename}, '
                             f'{counter.value}/{num_files_to_process}')


def clustalo_align_multiprocessing(fasta_to_align_folder,
                                   pool_threads=1,
                                   clustalo_threads=1,
                                   logger=None):
    """
    Generate alignments via function <clustalo_align> using multiprocessing.

    :param fasta_to_align_folder: path to folder containing input fasta alignment files
    :param int pool_threads: number of alignments to run concurrently
    :param int clustalo_threads: number of threads to use for each concurrent alignment
    :param logging.Logger logger: a logger object
    :return str output_folder: name of the output folder containing alignments
    """

    output_folder = f'11_pre_paralog_resolution_alignments'
    utils.createfolder(output_folder)

    logger.info(f'\n{"[INFO]:":10} Generating alignments for fasta files using Clustal Omega...')
    target_genes = [file for file in sorted(glob.glob(f'{fasta_to_align_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(clustalo_align,
                                      fasta_file,
                                      output_folder,
                                      counter,
                                      lock,
                                      num_files_to_process=len(target_genes),
                                      threads=clustalo_threads,
                                      logger=logger)
                          for fasta_file in target_genes]

        for future in as_completed(future_results):

            try:
                check = future.result()

            except Exception as error:
                print(f'Error raised: {error}')
                tb = traceback.format_exc()
                print(f'traceback is:\n{tb}')

        wait(future_results, return_when="ALL_COMPLETED")

    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                      utils.file_exists_and_not_empty(alignment)]

    logger.debug(f'{len(alignment_list)} alignments generated from {len(future_results)} fasta files...')

    return output_folder


def clustalo_align(fasta_file,
                   output_folder,
                   counter,
                   lock,
                   num_files_to_process,
                   threads=1,
                   logger=None):
    """
    Use clustal omega to align a fasta file of sequences, using the number of threads provided. Trims alignment with
    Trimal. Returns filename of the alignment produced.

    :param str fasta_file: path to a fasta file
    :param str output_folder: name of output folder for alignments
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int num_files_to_process: total number of fasta files for alignment
    :param int threads: number of threads to use for alignment program
    :param logging.Logger logger: a logger object
    :return str expected_alignment_file: filename of output alignment
    """

    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{fasta_file_basename}'

    try:
        assert utils.file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
        with lock:
            counter.value += 1

        return os.path.basename(expected_alignment_file)

    except AssertionError:
        command = f'clustalo --in {fasta_file} --out {expected_alignment_file} --threads {threads} --verbose'
        logger.info(f'{"[INFO]:":10} Performing Clustal alignment with command: {command}')

        try:

            result = subprocess.run(command,
                                    universal_newlines=True,
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.PIPE,
                                    check=True,
                                    shell=True)

            logger.debug(f'ClustalO check_returncode() is: {result.check_returncode()}')
            logger.debug(f'ClustalO stdout is: {result.stdout}')
            logger.debug(f'ClustalO stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'ClustalO FAILED. Output is: {exc}')
            logger.error(f'ClustalO stdout is: {exc.stdout}')
            logger.error(f'ClustalO stderr is: {exc.stderr}')

            raise ValueError('There was an issue running ClustalO. Check input files!')

        with lock:
            counter.value += 1
            logger.debug(f'Aligned file {fasta_file_basename}')

        return os.path.basename(expected_alignment_file)

    finally:
        with lock:
            sys.stderr.write(f'\r{"[INFO]:":10} Finished generating alignment for file {fasta_file_basename}, '
                             f'{counter.value}/{num_files_to_process}')


def fasttree_multiprocessing(alignments_folder,
                             pool=1,
                             threads=1,
                             bootstraps=False,
                             logger=None):
    """
    Generate FastTree trees using multiprocessing.

    :param str alignments_folder: name of folder containing alignments
    :param int pool: number of concurrent trees to run; default is 1
    :param int threads: number threads for each concurrent tree; default is 1
    :param bool bootstraps: if False, turn off support calculation
    :param logging.Logger logger: a logger object
    :return str output_folder: path to output folder containing tree files
    """

    output_folder = f'13_pre_paralog_resolution_trees'
    utils.createfolder(output_folder)

    logger.info(f'\n{"[INFO]:":10} Generating phylogenies from alignments using FastTreeMP...')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(fasttree,
                                      alignment,
                                      output_folder,
                                      threads,
                                      counter,
                                      lock,
                                      num_files_to_process=len(alignments),
                                      bootstraps=bootstraps,
                                      logger=logger)
                          for alignment in alignments]

        for future in as_completed(future_results):

            try:
                check = future.result()

            except Exception as error:
                print(f'Error raised: {error}')
                tb = traceback.format_exc()
                print(f'traceback is:\n{tb}')

        wait(future_results, return_when="ALL_COMPLETED")

    tree_list = [tree for tree in glob.glob(f'{output_folder}/*.treefile') if utils.file_exists_and_not_empty(tree)]

    logger.debug(f'{len(tree_list)} trees generated from {len(future_results)} fasta files...')

    return output_folder


def fasttree(alignment_file,
             output_folder,
             threads,
             counter,
             lock,
             num_files_to_process,
             bootstraps=False,
             logger=None):
    """
    Generate trees from alignments using FastTreeMP

    :param str alignment_file: path to alignment file
    :param str output_folder: name of output folder for tree
    :param int threads: number of threads to use for FastTreeMP
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int num_files_to_process: total number of alignment fasta files for tree generation
    :param bool bootstraps: if False, turn off support calculation
    :param logging.Logger logger: a logger object
    :return str expected_output_file: name of expected output tree file
    """

    alignment_file_basename = os.path.basename(alignment_file)
    expected_output_file = f'{output_folder}/{alignment_file_basename}.treefile'

    try:
        assert utils.file_exists_and_not_empty(expected_output_file)
        logger.debug(f'Output exists for {expected_output_file}, skipping...')
        with lock:
            counter.value += 1

        return os.path.basename(expected_output_file)

    except AssertionError:
        try:
            if bootstraps:
                fasttree_command = f'export OMP_NUM_THREADS={threads};' \
                                   f'FastTreeMP -gtr -nt < {alignment_file} > {expected_output_file}'
                result = subprocess.run(fasttree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
                logger.debug(f'FastTreeMP check_returncode() is: {result.check_returncode()}')
                logger.debug(f'FastTreeMP stdout is: {result.stdout}')
                logger.debug(f'FastTreeMP stderr is: {result.stderr}')

            else:
                fasttree_command = f'export OMP_NUM_THREADS={threads};' \
                                   f'FastTreeMP -gtr -nt -nosupport < {alignment_file} > {expected_output_file}'
                result = subprocess.run(fasttree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
                logger.debug(f'FastTreeMP check_returncode() is: {result.check_returncode()}')
                logger.debug(f'FastTreeMP stdout is: {result.stdout}')
                logger.debug(f'FastTreeMP stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'FastTreeMP FAILED. Output is: {exc}')
            logger.error(f'FastTreeMP stdout is: {exc.stdout}')
            logger.error(f'FastTreeMP stderr is: {exc.stderr}')

        with lock:
            counter.value += 1

        return os.path.basename(expected_output_file)

    finally:
        sys.stderr.write(f'\r{"[INFO]:":10} Finished generating tree {os.path.basename(expected_output_file)},'
                         f' {counter.value}/{num_files_to_process}')


def iqtree_multiprocessing(alignments_folder,
                           pool=1,
                           threads=1,
                           bootstraps=False,
                           logger=None):
    """
    Generate IQTree trees using multiprocessing.

    :param str alignments_folder: name of folder containing alignments
    :param int pool: number of concurrent trees to run; default is 1
    :param int threads: number threads for each concurrent tree; default is 1
    :param bool bootstraps: if False, turn off support calculation
    :param logging.Logger logger: a logger object
    :return str output_folder: path to output folder containing tree files
    """

    output_folder = f'13_pre_paralog_resolution_trees'
    utils.createfolder(output_folder)

    logger.info(f'\n{"[INFO]:":10} Generating phylogenies from alignments using IQTREE...')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*.fasta'))]

    with ProcessPoolExecutor(max_workers=pool) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(iqtree,
                                      alignment,
                                      output_folder,
                                      threads,
                                      counter,
                                      lock,
                                      num_files_to_process=len(alignments),
                                      bootstraps=bootstraps,
                                      logger=logger)
                          for alignment in alignments]

        for future in as_completed(future_results):

            try:
                check = future.result()

            except Exception as error:
                print(f'Error raised: {error}')
                tb = traceback.format_exc()
                print(f'traceback is:\n{tb}')

        wait(future_results, return_when="ALL_COMPLETED")

    tree_list = [tree for tree in glob.glob(f'{output_folder}/*.treefile') if utils.file_exists_and_not_empty(tree)]

    logger.debug(f'\n{len(tree_list)} trees generated from {len(future_results)} fasta files...')

    return output_folder


def iqtree(alignment_file,
           output_folder,
           threads,
           counter,
           lock,
           num_files_to_process,
           bootstraps=False,
           logger=None):
    """
    Generate trees from alignments using IQTREE

    :param str alignment_file: path to alignment file
    :param str output_folder: name of output folder for tree
    :param int threads: number of threads to use for FastTreeMP
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int num_files_to_process: total number of alignment fasta files for tree generation
    :param bool bootstraps: if False, turn off support calculation
    :param logging.Logger logger: a logger object
    :return str expected_output_file: name of expected output tree file
    """

    alignment_file_basename = os.path.basename(alignment_file)
    expected_output_file = f'{output_folder}/{alignment_file_basename}.treefile'

    try:
        assert utils.file_exists_and_not_empty(expected_output_file)
        logger.debug(f'Output exists for {expected_output_file}, skipping...')
        with lock:
            counter.value += 1

        return os.path.basename(expected_output_file)

    except AssertionError:
        try:
            if bootstraps:
                iqtree_command = f'iqtree -redo -pre {output_folder}/{alignment_file_basename} -s {alignment_file} ' \
                                 f'-m GTR+G -bb 1000 -bnni -nt {str(threads)} -quiet'
                result = subprocess.run(iqtree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
                logger.debug(f'IQTREE check_returncode() is: {result.check_returncode()}')
                logger.debug(f'IQTREE stdout is: {result.stdout}')
                logger.debug(f'IQTREE stderr is: {result.stderr}')

            else:
                iqtree_command = f'iqtree -redo -pre {output_folder}/{alignment_file_basename} -s {alignment_file} ' \
                                 f'-m GTR+G -nt {str(threads)} -quiet'
                result = subprocess.run(iqtree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True, check=True)
                logger.debug(f'IQTREE check_returncode() is: {result.check_returncode()}')
                logger.debug(f'IQTREE stdout is: {result.stdout}')
                logger.debug(f'IQTREE stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'IQTREE FAILED. Output is: {exc}')
            logger.error(f'IQTREE stdout is: {exc.stdout}')
            logger.error(f'IQTREE stderr is: {exc.stderr}')

        with lock:
            counter.value += 1

        return os.path.basename(expected_output_file)

    finally:
        sys.stderr.write(f'\r{"[INFO]:":10} Finished generating tree {os.path.basename(expected_output_file)},'
                         f' {counter.value}/{num_files_to_process}')


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

    logger.debug(f'{"[INFO]:":10} Module align_selected_and_tree was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.debug(f'{fill}\n')
    logger.debug(args)

    logger.info('')
    logger.info(f'{"[INFO]:":10} ======> ALIGNMENT AND TREE FROM SELECTED SEQUENCES <======\n')

    # Checking input directories and files:
    selected_alignments_directory = '09_sequences_from_qc_trees'
    selected_alignments_suffix = 'selected.fasta'
    qc_alignments_directory = args.qc_alignment_directory  # i.e. with all seqs, before tree quality control
    qc_alignments_suffix = '.fasta'

    directory_suffix_dict = {selected_alignments_directory: selected_alignments_suffix,
                             qc_alignments_directory: qc_alignments_suffix}
    file_list = []

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Get a list of internal and external taxon names from log file "outgroup_taxon_list.tsv":
    try:
        internal_outgroup_taxa = []
        external_outgroup_taxa = []

        with open(f'00_logs_and_reports/reports/outgroup_taxon_list.tsv', 'r') as outgroup_taxon_handle:
            lines = outgroup_taxon_handle.readlines()
            for line in lines:
                outgroup_type, taxon = line.split()
                if outgroup_type == 'INTERNAL_OUTGROUP':
                    internal_outgroup_taxa.append(taxon)
                if outgroup_type == 'EXTERNAL_OUTGROUP':
                    external_outgroup_taxa.append(taxon)

        fill = utils.fill_forward_slash(f'{"[INFO]:":10} External outgroup taxa: {external_outgroup_taxa}',
                                        width=90, subsequent_indent=' ' * 11, break_on_forward_slash=True)
        logger.info(f'{fill}')

        fill = utils.fill_forward_slash(f'{"[INFO]:":10} Internal outgroup taxa: {internal_outgroup_taxa}',
                                        width=90, subsequent_indent=' ' * 11, break_on_forward_slash=True)
        logger.info(f'{fill}')

    except FileNotFoundError:
        logger.error(f'{"[ERROR]:":10} No outgroup taxon list found at '    
                     f'00_logs_and_reports/reports/outgroup_taxon_list.tsv')
        sys.exit()

    # Add outgroup sequences, both internal (if removed by the tree QC steps) and external (if a fasta file of
    # external outgroup sequences is provided):

    outgroups_added_folder = add_outgroup_seqs(args.qc_alignment_directory,
                                               selected_alignments_directory,
                                               internal_outgroup_taxa,
                                               external_outgroup_taxa,
                                               logger=logger)

    if not args.use_clustal:
        logger.debug(f'Running without --use_clustal option - aligning with MAFFT only')

        alignments_output_folder = mafft_align_multiprocessing(
            outgroups_added_folder,
            algorithm=args.mafft_algorithm,
            pool_threads=args.pool,
            mafft_threads=args.threads,
            logger=logger)

        # Perform optional trimming with TrimAl:
        trimmed_output_folder = f'12_pre_paralog_resolution_alignments_trimmed'
        if not args.no_trimming:
            alignments_output_folder = run_trimal(alignments_output_folder,
                                                  trimmed_output_folder,
                                                  args_for_trimal_options=args,
                                                  logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping trimming step...')

        # Generate trees:
        if args.use_fasttree:
            trees_folder = fasttree_multiprocessing(alignments_output_folder,
                                                    pool=args.pool,
                                                    threads=args.threads,
                                                    bootstraps=args.generate_bootstraps,
                                                    logger=logger)

            utils.resolve_polytomies(trees_folder, logger=logger)

        else:
            trees_folder = iqtree_multiprocessing(alignments_output_folder,
                                                  pool=args.pool,
                                                  threads=args.threads,
                                                  bootstraps=args.generate_bootstraps,
                                                  logger=logger)

            utils.resolve_polytomies(trees_folder, logger=logger)

    elif args.use_clustal:  # Re-align with Clustal Omega only - no need for MAFFT.
        logger.debug(f'Running with --use-clustal option - aligning with Clustal Omega only')

        alignments_output_folder = clustalo_align_multiprocessing(
            outgroups_added_folder,
            pool_threads=args.pool,
            clustalo_threads=args.threads,
            logger=logger)

        # Perform optional trimming with TrimAl:
        trimmed_output_folder = f'12_pre_paralog_resolution_alignments_trimmed'
        if not args.no_trimming:
            alignments_output_folder = run_trimal(alignments_output_folder,
                                                  trimmed_output_folder,
                                                  args_for_trimal_options=args,
                                                  logger=logger)
        else:
            logger.info(f'\n{"[INFO]:":10} Skipping trimming step...')

        # Generate trees:
        if args.use_fasttree:
            trees_folder = fasttree_multiprocessing(alignments_output_folder,
                                                    pool=args.pool,
                                                    threads=args.threads,
                                                    bootstraps=args.generate_bootstraps,
                                                    logger=logger)

            utils.resolve_polytomies(trees_folder, logger=logger)

        else:
            trees_folder = iqtree_multiprocessing(alignments_output_folder,
                                                  pool=args.pool,
                                                  threads=args.threads,
                                                  bootstraps=args.generate_bootstraps,
                                                  logger=logger)

            utils.resolve_polytomies(trees_folder, logger=logger)

    logger.info(f'{"[INFO]:":10} Finished aligning selected fasta sequences and generating trees.')
