#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au https://github.com/chrisjackson-pellicle

"""
Input is fasta files that have undergone QC processes; this QC may have removed internal outgroups (if specified).
This script:

    - Checks for any paralogs in internal outgroups (if the latter are specified), and selects a single
    representative sequence for each taxon  - this will be the sequence with the highest distance from a sample of up
    to 10 ingroup sequences.
    - Aligns the sequences using either MAFFT or MUSCLE.
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
from Bio.Align.Applications import MafftCommandline, ClustalOmegaCommandline, MuscleCommandline
from concurrent.futures.process import ProcessPoolExecutor
from multiprocessing import Manager
from concurrent.futures import wait

from yang_and_smith import utils


def add_outgroup_seqs(hmmcleaned_alignment_directory,
                      selected_alignment_directory,
                      list_of_internal_outgroups,
                      file_of_external_outgroups,
                      list_of_external_outgroups=None,
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

    :param str hmmcleaned_alignment_directory:
    :param str selected_alignment_directory:
    :param list list_of_internal_outgroups:
    :param str file_of_external_outgroups:
    :param list_of_external_outgroups:
    :param logging.Logger logger: a logger object
    :return str output_folder: name of output folder containing fasta with outgroups added
    """

    logger.debug(f'list_of_internal_outgroups: {list_of_internal_outgroups}')
    logger.debug(f'list_of_external_outgroups_to_select: {list_of_external_outgroups}')

    if not file_of_external_outgroups and not list_of_internal_outgroups:
        logger.warning(f'{"[WARNING]:":10} No external or internal outgroups supplied!')

    input_folder_basename = os.path.basename(selected_alignment_directory)
    output_folder = f'13_{input_folder_basename.lstrip("12_")}_selected_alignments_outgroups_added'
    utils.createfolder(output_folder)  # for the outgroups added fasta files

    # Read in original paralog fasta files, and create a dictionary of gene_id:list_of_seq_names for taxa in
    # list_of_internal_outgroups:
    internal_outgroup_dict = defaultdict(lambda: defaultdict(list))
    all_paralog_taxon_names = set()

    for original_alignment in glob.glob(f'{hmmcleaned_alignment_directory}/*.hmm.trimmed.fasta'):
        gene_id = os.path.basename(original_alignment).split('.')[0]  # get prefix e.g. '4471'
        seqs = SeqIO.parse(original_alignment, 'fasta')

        # Populate the internal_outgroup_dict (potentially more than one sequence per taxon):
        for seq in seqs:
            seq_name_prefix = seq.name.split('.')[0]  # assumes there are no other dots ('.') in the sequence name
            all_paralog_taxon_names.add(seq_name_prefix)
            if list_of_internal_outgroups and seq_name_prefix in list_of_internal_outgroups:
                internal_outgroup_dict[gene_id][seq_name_prefix].append(seq)

        if list_of_internal_outgroups:
            alignment = AlignIO.read(original_alignment, 'fasta')

            # Create an MultipleSeqAlignment object for ingroup sequences only, and take the first 10 sequences from
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

    # Read in external outgroups file if present, and create a dictionary of gene_id:list_of_seq_names, either for all
    # seqs if no external outgroup taxa specified, or for specified taxa only:
    external_outgroup_dict = defaultdict(list)
    if file_of_external_outgroups:  # dict not populated if no outgroups file provided
        seqs = SeqIO.parse(file_of_external_outgroups, 'fasta')
        for seq in seqs:
            gene_id = seq.name.split('-')[-1]  # get gene id e.g. '4471'
            taxon = '-'.join(seq.name.split('-')[:-1])  # e.g. 'AMBTR'
            if list_of_external_outgroups:
                if taxon in list_of_external_outgroups:
                    seq.name = taxon  # i.e. we don't want the suffix e.g. '4471' present in the fasta file
                    seq.id = taxon
                    external_outgroup_dict[gene_id].append(seq)
            else:
                seq.name = taxon  # i.e. we don't want the suffix e.g. '4471' present in the fasta file
                seq.id = taxon
                external_outgroup_dict[gene_id].append(seq)

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
        logger.info(f'ingroup_taxon_names: {ingroup_taxon_names}')
    else:
        ingroup_taxon_names = [name for name in all_paralog_taxon_names]
        logger.debug(f'ingroup_taxon_names: {ingroup_taxon_names}')
    with open(f'in_and_outgroups_list_{os.path.basename(selected_alignment_directory)}.txt', 'w') as group_list:
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
            logger.info(f'Taxon {taxon_id} is an internal outgroup, and has more than one sequence for gene '
                        f'{gene_id}. Only the sequence most divergent from the ingroup taxa will be retained.')

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
                if closest_sequence_distance > seq_to_keep_distance:
                    logger.debug(f'{closest_sequence_distance} is larger than {seq_to_keep_distance}')
                    seq_to_keep_distance = closest_sequence_distance
                    seq_to_keep = seq
            logger.debug(f'Keeping sequence {seq_to_keep.name} with distance {seq_to_keep_distance}')
            internal_outgroup_dict_copy[gene_id][taxon_id] = [seq_to_keep]

    return internal_outgroup_dict_copy


def mafft_or_muscle_align_multiprocessing(fasta_to_align_folder,
                                          algorithm='linsi',
                                          pool_threads=1,
                                          mafft_threads=1,
                                          no_stitched_contigs=False,
                                          use_muscle=False,
                                          logger=None):

    """
    Generate alignments via function <mafft_or_muscle_align> using multiprocessing.

    :param str fasta_to_align_folder: path to folder containing input fasta files with outgroups added
    :param str algorithm: algorithm to use for mafft alignment; default is 'auto'
    :param int pool_threads: number of alignments to run concurrently
    :param int mafft_threads: number of threads to use for each concurrent alignment
    :param bool no_stitched_contigs: if True, realign with Clustal Omega
    :param bool use_muscle: if True, use muscle instead of mafft for alignments
    :param logging.Logger logger: a logger object
    :return str output_folder: name of the output folder containing alignments
    """

    input_folder_basename = os.path.basename(fasta_to_align_folder)
    output_folder = f'14_{input_folder_basename.lstrip("13_")}_alignments'
    utils.createfolder(output_folder)

    if use_muscle:
        logger.info(f'{"[INFO]:":10} Generating alignments for fasta files in folder {input_folder_basename} using '
                    f'MUSCLE...')
    else:
        logger.info(f'{"[INFO]:":10} Generating alignments for fasta files in folder {input_folder_basename} using '
                    f'MAFFT...')

    # Filter out any input files with fewer than four sequences:
    target_genes = []
    for fasta_file in sorted(glob.glob(f'{fasta_to_align_folder}/*.outgroup_added.fasta')):
        with open(fasta_file, 'r') as input_fasta_handle:
            seqs = list(SeqIO.parse(input_fasta_handle, 'fasta'))
            if len(seqs) < 4:
                logger.info(f'{"[INFO]:":10} Skipping file {fasta_file} as it contains fewer than 4 sequences!')
                continue
            else:
                target_genes.append(fasta_file)

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(mafft_or_muscle_align,
                                      fasta_file,
                                      algorithm,
                                      output_folder,
                                      counter,
                                      lock,
                                      num_files_to_process=len(target_genes),
                                      threads=mafft_threads,
                                      no_stitched_contigs=no_stitched_contigs,
                                      use_muscle=use_muscle,
                                      logger=logger)

                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(utils.done_callback)
        wait(future_results, return_when="ALL_COMPLETED")

    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                      utils.file_exists_and_not_empty(alignment)]

    logger.debug(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')

    return output_folder


def mafft_or_muscle_align(fasta_file,
                          algorithm,
                          output_folder,
                          counter,
                          lock,
                          num_files_to_process,
                          threads=1,
                          no_stitched_contigs=False,
                          use_muscle=False,
                          logger=None):
    """
    Use mafft or muscle to align a fasta file of sequences, using the algorithm (if mafft) and number of threads
    provided. Trims alignment with Trimal if no_stitched_contigs=False. Returns filename of the output alignment.

    :param str fasta_file: path to a fasta file
    :param str algorithm: algorithm to use for mafft alignment; default is 'auto'
    :param str output_folder: name of output folder for alignments
    :param multiprocessing.managers.ValueProxy counter: shared counter for fasta files processed
    :param multiprocessing.managers.AcquirerProxy lock: lock for ordered logging of info messages
    :param int num_files_to_process: total number of fasta files for alignment
    :param int threads: number of threads to use for alignment program
    :param bool no_stitched_contigs: if True, realign with Clustal Omega
    :param bool use_muscle: if True, use muscle instead of mafft for alignments
    :param logging.Logger logger: a logger object
    :return str expected_alignment_file/expected_alignment_file_trimmed: filename of output alignment
    """

    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'
    expected_alignment_file_trimmed = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)

    try:
        if not no_stitched_contigs:  # i.e. stitched contig _were_ produced
            assert utils.file_exists_and_not_empty(expected_alignment_file_trimmed)
            logger.debug(f'Trimmed alignment exists for {fasta_file_basename}, skipping...')
            with lock:
                counter.value += 1

            return os.path.basename(expected_alignment_file_trimmed)
        else:
            assert utils.file_exists_and_not_empty(expected_alignment_file)
            logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
            with lock:
                counter.value += 1

            return os.path.basename(expected_alignment_file)

    except AssertionError:
        if use_muscle:

            logger.info(f'{"[INFO]:":10} Alignment will be performed using MUSCLE rather than MAFFT!')
            muscle_cline = MuscleCommandline(input=fasta_file, out=expected_alignment_file)
            stdout, stderr = muscle_cline()

            logger.debug(f'stdout is: {stdout}')
            logger.debug(f'stderr is: {stderr}')
        else:
            if algorithm == 'auto':
                mafft_cline = (MafftCommandline(auto='true', thread=threads, input=fasta_file))
            else:
                mafft_cline = (MafftCommandline(algorithm, thread=threads, input=fasta_file))

            logger.info(f'{"[INFO]:":10} Performing MAFFT alignment with command: {mafft_cline}')
            stdout, stderr = mafft_cline()

            with open(expected_alignment_file, 'w') as alignment_file:
                alignment_file.write(stdout)

        if not no_stitched_contigs:  # only trim if stitched contigs (and hence won't be realigned)
            trimmed_alignment = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)
            try:
                result = subprocess.run(['trimal', '-in', expected_alignment_file, '-out', trimmed_alignment,
                                         '-gapthreshold', '0.12', '-terminalonly', '-gw', '1'],
                                        universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        check=True)
                logger.debug(f'trimal check_returncode() is: {result.check_returncode()}')
                logger.debug(f'trimal stdout is: {result.stdout}')
                logger.debug(f'trimal stderr is: {result.stderr}')

            except subprocess.CalledProcessError as exc:
                logger.error(f'trimal FAILED. Output is: {exc}')
                logger.error(f'trimal stdout is: {exc.stdout}')
                logger.error(f'trimal stderr is: {exc.stderr}')
                raise ValueError('There was an issue running trimal. Check input files!')

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

    input_folder_basename = os.path.basename(fasta_to_align_folder)
    output_folder = f'{input_folder_basename}_clustal'
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
        for future in future_results:
            future.add_done_callback(utils.done_callback)
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
        clustalomega_cline = ClustalOmegaCommandline(infile=fasta_file, outfile=expected_alignment_file,
                                                     verbose=True, auto=True, threads=threads)
        logger.info(f'{"[INFO]:":10} Performing Clustal alignment with command: {clustalomega_cline}')

        stdout, stderr = clustalomega_cline()

        logger.debug(f'stdout is: {stdout}')
        logger.debug(f'stderr is: {stderr}')

        trimmed_alignment = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)

        try:

            result = subprocess.run(['trimal', '-in', expected_alignment_file, '-out', trimmed_alignment,
                                     '-gapthreshold', '0.12', '-terminalonly', '-gw', '1'],
                                    universal_newlines=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                    check=True)
            logger.debug(f'trimal check_returncode() is: {result.check_returncode()}')
            logger.debug(f'trimal stdout is: {result.stdout}')
            logger.debug(f'trimal stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'trimal FAILED. Output is: {exc}')
            logger.error(f'trimal stdout is: {exc.stdout}')
            logger.error(f'trimal stderr is: {exc.stderr}')
            raise ValueError('There was an issue running trimal. Check input files!')

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

    input_folder_basename = os.path.basename(alignments_folder)
    output_folder = f'15_{input_folder_basename.lstrip("14_")}_tree_files'
    utils.createfolder(output_folder)

    logger.info(f'\n{"[INFO]:":10} Generating phylogenies from alignments using FastTreeMP...')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*trimmed.fasta'))]

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
        for future in future_results:
            future.add_done_callback(utils.done_callback)

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

    input_folder_basename = os.path.basename(alignments_folder)
    output_folder = f'15_{input_folder_basename.lstrip("14_")}_tree_files'
    utils.createfolder(output_folder)

    logger.info(f'\n{"[INFO]:":10} Generating phylogenies from alignments using IQTREE...')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*.trimmed.fasta'))]

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
        for future in future_results:
            future.add_done_callback(utils.done_callback)

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


def main(args):
    """
    Entry point for the resolve_paralogs.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = utils.setup_logger(__name__, '00_logs_and_reports_resolve_paralogs/logs/09_align_selected_and_tree')

    # check for external dependencies:
    if utils.check_dependencies(logger=logger):
        logger.info(f'{"[INFO]:":10} All external dependencies found!')
    else:
        logger.error(f'{"[ERROR]:":10} One or more dependencies not found!')
        sys.exit(1)

    logger.info(f'{"[INFO]:":10} Subcommand align_selected_and_tree was called with these arguments:')
    fill = textwrap.fill(' '.join(sys.argv[1:]), width=90, initial_indent=' ' * 11, subsequent_indent=' ' * 11,
                         break_on_hyphens=False)
    logger.info(f'{fill}\n')
    logger.debug(args)

    # Checking input directories and files:
    directory_suffix_dict = {args.selected_alignment_directory: '.selected.fasta',
                             args.hmmcleaned_alignment_directory: 'aln.hmm.trimmed.fasta'}
    file_list = []

    if args.external_outgroups_file:
        file_list.append(args.external_outgroups_file)

    utils.check_inputs(directory_suffix_dict,
                       file_list,
                       logger=logger)

    # Add outgroup sequences, both internal (if removed by the tree QC steps) and external (if a fasta file of
    # external outgroup sequences is provided).
    outgroups_added_folder = add_outgroup_seqs(args.hmmcleaned_alignment_directory,
                                               args.selected_alignment_directory,
                                               args.internal_outgroups,
                                               args.external_outgroups_file,
                                               list_of_external_outgroups=args.external_outgroups,
                                               logger=logger)

    if not args.no_stitched_contigs:  # i.e. if it's a standard run with stitched contigs produced.
        logger.debug(f'Running without no_stitched_contigs option - aligning with mafft or muscle only')

        alignments_output_folder = mafft_or_muscle_align_multiprocessing(
            outgroups_added_folder,
            algorithm=args.mafft_algorithm,
            pool_threads=args.pool,
            mafft_threads=args.threads,
            no_stitched_contigs=args.no_stitched_contigs,
            use_muscle=args.use_muscle,
            logger=logger)

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

    elif args.no_stitched_contigs:  # re-align with Clustal Omega.
        logger.debug(f'Running with no_stitched_contigs option - realigning with clustal omega')

        alignments_output_folder = mafft_or_muscle_align_multiprocessing(
            outgroups_added_folder,
            algorithm=args.mafft_algorithm,
            pool_threads=args.pool,
            mafft_threads=args.threads,
            no_stitched_contigs=args.no_stitched_contigs,
            use_muscle=args.use_muscle,
            logger=logger)

        clustal_alignment_output_folder = clustalo_align_multiprocessing(
            alignments_output_folder,
            pool_threads=args.pool,
            clustalo_threads=args.threads,
            logger=logger)

        # Generate trees:
        if args.use_fasttree:
            trees_folder = fasttree_multiprocessing(clustal_alignment_output_folder,
                                                    pool=args.pool,
                                                    threads=args.threads,
                                                    bootstraps=args.generate_bootstraps,
                                                    logger=logger)

            utils.resolve_polytomies(trees_folder, logger=logger)

        else:
            trees_folder = iqtree_multiprocessing(clustal_alignment_output_folder,
                                                  pool=args.pool,
                                                  threads=args.threads,
                                                  bootstraps=args.generate_bootstraps,
                                                  logger=logger)

            utils.resolve_polytomies(trees_folder, logger=logger)

    logger.info(f'{"[INFO]:":10} Finished aligning selected fasta sequences and generating trees.')
