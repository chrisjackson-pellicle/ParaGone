#!/usr/bin/env python

# Author: Chris Jackson chris.jackson@rbg.vic.gov.au

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




def mafft_align(fasta_file, algorithm, output_folder, counter, lock, num_files_to_process, threads=2,
                no_supercontigs=False, use_muscle=False):
    """
    Uses mafft to align a fasta file of sequences, using the algorithm and number of threads provided. Returns filename
    of the alignment produced.
    """

    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{re.sub(".fasta", ".aln.fasta", fasta_file_basename)}'
    expected_alignment_file_trimmed = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)
    # print(expected_alignment_file_trimmed)

    try:
        if not no_supercontigs:
            assert file_exists_and_not_empty(expected_alignment_file_trimmed)
            logger.debug(f'Trimmed alignment exists for {fasta_file_basename}, skipping...')
            with lock:
                counter.value += 1
            return os.path.basename(expected_alignment_file_trimmed)
        else:
            assert file_exists_and_not_empty(expected_alignment_file)
            logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
            with lock:
                counter.value += 1
            return os.path.basename(expected_alignment_file)
    except AssertionError:
        if use_muscle:
            # print(fasta_file_basename)
            logger.info('Alignment will be performed using MUSCLE rather than MAFFT!')
            muscle_cline = MuscleCommandline(input=fasta_file, out=expected_alignment_file)
            stdout, stderr = muscle_cline()
        else:
            if algorithm == 'auto':
                mafft_cline = (MafftCommandline(auto='true', thread=threads, input=fasta_file))
            else:
                mafft_cline = (MafftCommandline(algorithm, thread=threads, input=fasta_file))
            logger.info(f'Performing MAFFT alignment with command" {mafft_cline}')
            stdout, stderr = mafft_cline()
            with open(expected_alignment_file, 'w') as alignment_file:
                alignment_file.write(stdout)

        if not no_supercontigs:
            trimmed_alignment = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)
            # run_trim = subprocess.run(['/Users/chrisjackson/miniconda3/bin/trimal', '-in', expected_alignment_file, '-out', trimmed_alignment,
            #                            '-gapthreshold', '0.12', '-terminalonly', '-gw', '1'], check=True)
            run_trim = subprocess.run(['trimal', '-in', expected_alignment_file, '-out', trimmed_alignment,
                                       '-gapthreshold', '0.12', '-terminalonly', '-gw', '1'], check=True)
        with lock:
            counter.value += 1
        logger.debug(f'Aligned file {fasta_file_basename}')
        return os.path.basename(expected_alignment_file)
    finally:
        print(f'\rFinished generating alignment for file {fasta_file_basename}, '
              f'{counter.value}/{num_files_to_process}', end='')


def mafft_align_multiprocessing(fasta_to_align_folder, algorithm='linsi', pool_threads=1,
                                mafft_threads=2, no_supercontigs=False, use_muscle=False):
    """
    Generate alignments via function <align_targets> using multiprocessing.
    """

    input_folder_basename = os.path.basename(fasta_to_align_folder)
    output_folder = f'{input_folder_basename}_alignments'
    # print(f'output_folder from mafft_align_multiprocessing: {output_folder}')
    createfolder(output_folder)

    logger.info(f'Generating alignments for fasta files in folder {fasta_to_align_folder}...\n')
    # target_genes = [file for file in sorted(glob.glob(f'{fasta_to_align_folder}/*.outgroup_added.fasta'))]
    # print(target_genes)
    target_genes = []
    for fasta_file in sorted(glob.glob(f'{fasta_to_align_folder}/*.outgroup_added.fasta')):
        with open(fasta_file, 'r') as input_fasta_handle:
            seqs = list(SeqIO.parse(input_fasta_handle, 'fasta'))
            if len(seqs) < 4:
                logger.warning(f'Skipping file {fasta_file} as it contains fewer than 4 sequences!')
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
                                      no_supercontigs=no_supercontigs,
                                      use_muscle=use_muscle)

                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                      file_exists_and_not_empty(alignment)]
    logger.info(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')

    return output_folder


def clustalo_align(fasta_file, output_folder, counter, lock, num_files_to_process, threads=2):
    """
    Uses clustal omega to align a fasta file of sequences, using the algorithm and number of threads provided. Returns
    filename
    of the alignment produced.
    """

    createfolder(output_folder)
    fasta_file_basename = os.path.basename(fasta_file)
    expected_alignment_file = f'{output_folder}/{fasta_file_basename}'

    try:
        assert file_exists_and_not_empty(expected_alignment_file)
        logger.debug(f'Alignment exists for {fasta_file_basename}, skipping...')
    except AssertionError:
        # clustalomega_cline = ClustalOmegaCommandline(infile=fasta_file, outfile=expected_alignment_file,
        #                                              verbose=True, auto=True, threads=threads, pileup=True)
        clustalomega_cline = ClustalOmegaCommandline(infile=fasta_file, outfile=expected_alignment_file,
                                                     verbose=True, auto=True, threads=threads)
        clustalomega_cline()

        trimmed_alignment = re.sub('.aln.fasta', '.aln.trimmed.fasta', expected_alignment_file)
        run_trim = subprocess.run(['trimal', '-in', expected_alignment_file, '-out', trimmed_alignment, '-gapthreshold',
                                   '0.12', '-terminalonly', '-gw', '1'], check=True)
    finally:
        with lock:
            counter.value += 1
            print(f'\rFinished generating alignment for file {fasta_file_basename}, '
                  f'{counter.value}/{num_files_to_process}', end='')
        return os.path.basename(expected_alignment_file)


def clustalo_align_multiprocessing(fasta_to_align_folder, pool_threads=1, clustalo_threads=2):
    """
    Generate alignments via function <clustalo_align> using multiprocessing.
    """

    input_folder_basename = os.path.basename(fasta_to_align_folder)
    output_folder = f'{input_folder_basename}_clustal'
    # print(f'output_folder from clustal_align_multiprocessing: {output_folder}')
    createfolder(output_folder)

    logger.debug('Generating alignments for fasta files using clustal omega...\n')
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
                                      threads=clustalo_threads)
                          for fasta_file in target_genes]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    alignment_list = [alignment for alignment in glob.glob(f'{output_folder}/*.aln.fasta') if
                      file_exists_and_not_empty(alignment)]
    logger.info(f'\n{len(alignment_list)} alignments generated from {len(future_results)} fasta files...\n')

    return output_folder


def iqtree(alignment_file, output_folder, iqtree_threads, counter, lock, num_files_to_process, bootstraps=False):
    """
    Generate trees from alignments using iqtree
    """
    alignment_file_basename = os.path.basename(alignment_file)
    expected_output_file = f'{output_folder}/{alignment_file_basename}.treefile'
    print(alignment_file)

    try:
        assert file_exists_and_not_empty(expected_output_file)
        logger.debug(f'Output exists for {expected_output_file}, skipping...')
        with lock:
            counter.value += 1
        return os.path.basename(expected_output_file)
    except AssertionError:
        try:
            if bootstraps:
                check_iqtree = subprocess.run(['iqtree', '-redo', '-pre',
                                               f'{output_folder}/{alignment_file_basename}',
                                               '-s', alignment_file, '-m', 'GTR+G', '-bb', '1000', '-bnni', '-nt',
                                               str(iqtree_threads), '-quiet'], check=True)
            else:
                check_iqtree = subprocess.run(['iqtree', '-redo', '-pre',
                                               f'{output_folder}/{alignment_file_basename}',
                                               '-s', alignment_file, '-m', 'GTR+G', '-nt',
                                               str(iqtree_threads), '-quiet'], check=True)
                # print(check_iqtree)
        except:
            logger.info(f'No tree produced for {alignment_file}- fewer than 3 sequences in alignment?')
        with lock:
            counter.value += 1
        return os.path.basename(expected_output_file)
    finally:
        print(f'\rFinished generating output {os.path.basename(expected_output_file)}, {counter.value}/{num_files_to_process}', end='')


def iqtree_multiprocessing(alignments_folder, pool_threads=1, iqtree_threads=2, bootstraps=False):
    """
    Generate iqtree trees using multiprocessing.
    """

    input_folder_basename = os.path.basename(alignments_folder)
    output_folder = f'{input_folder_basename}_tree_files'
    createfolder(output_folder)

    logger.info('Generating trees from alignments...\n')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*.trimmed.fasta'))]
    # print(alignments)

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(iqtree, alignment, output_folder, iqtree_threads, counter, lock,
                                      num_files_to_process=len(alignments), bootstraps=bootstraps)
                          for alignment in alignments]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    tree_list = [tree for tree in glob.glob(f'{output_folder}/*.treefile') if file_exists_and_not_empty(tree)]
    logger.info(f'\n{len(tree_list)} trees generated from {len(future_results)} fasta files...\n')


def fasttree_multiprocessing(alignments_folder, pool_threads=1, bootstraps=False):
    """
    Generate FastTree trees using multiprocessing.
    """

    input_folder_basename = os.path.basename(alignments_folder)
    output_folder = f'{input_folder_basename}_tree_files'
    createfolder(output_folder)

    logger.info('Generating FastTree trees from alignments...\n')
    alignments = [file for file in sorted(glob.glob(f'{alignments_folder}/*trimmed.fasta'))]
    alignments.extend([file for file in sorted(glob.glob(f'{alignments_folder}/*trimmed_hmm.fasta'))])
    # print(alignments)

    with ProcessPoolExecutor(max_workers=pool_threads) as pool:
        manager = Manager()
        lock = manager.Lock()
        counter = manager.Value('i', 0)
        future_results = [pool.submit(fasttree, alignment, output_folder, counter, lock,
                                      num_files_to_process=len(alignments), bootstraps=bootstraps)
                          for alignment in alignments]
        for future in future_results:
            future.add_done_callback(done_callback)
        wait(future_results, return_when="ALL_COMPLETED")
    tree_list = [tree for tree in glob.glob(f'{output_folder}/*.treefile') if file_exists_and_not_empty(tree)]
    logger.info(f'\n{len(tree_list)} trees generated from {len(future_results)} fasta files...\n')


def fasttree(alignment_file, output_folder, counter, lock, num_files_to_process, bootstraps=False):
    """
    Generate trees from alignments using FastTree
    """
    alignment_file_basename = os.path.basename(alignment_file)
    expected_output_file = f'{output_folder}/{alignment_file_basename}.treefile'

    try:
        assert file_exists_and_not_empty(expected_output_file)
        logger.debug(f'Output exists for {expected_output_file}, skipping...')
        with lock:
            counter.value += 1
        return os.path.basename(expected_output_file)
    except AssertionError:
        try:
            if bootstraps:
                fasttree_command = f'FastTreeMP -gtr -nt < {alignment_file} > {expected_output_file}'
                result = subprocess.run(fasttree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
                logger.debug(f'FastTreeMP check_returncode() is: {result.check_returncode()}')
                logger.debug(f'FastTreeMP stdout is: {result.stdout}')
                logger.debug(f'FastTreeMP stderr is: {result.stderr}')

            else:
                fasttree_command = f'FastTreeMP -gtr -nt -nosupport < {alignment_file} > {expected_output_file}'
                result = subprocess.run(fasttree_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                        universal_newlines=True)
                logger.debug(f'FastTreeMP check_returncode() is: {result.check_returncode()}')
                logger.debug(f'FastTreeMP stdout is: {result.stdout}')
                logger.debug(f'FastTreeMP stderr is: {result.stderr}')

        except subprocess.CalledProcessError as exc:
            logger.error(f'FastTreeMP FAILED. Output is: {exc}')
            logger.error(f'FastTreeMP stdout is: {exc.stdout}')
            logger.error(f'FastTreeMP stderr is: {exc.stderr}')

        # except:
        #     logger.info(f'\nNo tree produced for {alignment_file}- fewer than 3 sequences in alignment?\n')
        with lock:
            counter.value += 1
        return os.path.basename(expected_output_file)
    finally:
        print(f'\rFinished generating output {os.path.basename(expected_output_file)}, {counter.value}'
              f'/{num_files_to_process}', end='')


def filter_internal_outgroups(internal_outgroup_dict, alignment_ingroup_only, gene_id):
    """
    Filters internal outgroups with paralogs to retain the single sequence most divergent from the ingroup
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


def add_outgroup_seqs(original_paralog_gene_fasta_directory, folder_of_qc_paralog_files, list_of_internal_outgroups,
                      file_of_external_outgroups, list_of_external_outgroups=None):
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

    """

    logger.info(f'list_of_internal_outgroups: {list_of_internal_outgroups}')
    logger.info(f'list_of_external_outgroups: {list_of_external_outgroups}')

    input_folder_basename = os.path.basename(folder_of_qc_paralog_files)
    output_folder = f'{input_folder_basename}_outgroups_added'
    createfolder(output_folder)  # for the outgroups added fasta files

    # Read in original paralog fasta files, and create a dictionary of gene_id:list_of_seq_names for taxa in
    # list_of_internal_outgroups:
    internal_outgroup_dict = defaultdict(lambda: defaultdict(list))
    all_paralog_taxon_names = set()

    for fasta in glob.glob(f'{original_paralog_gene_fasta_directory}/*.hmm.trimmed.fasta'):
        gene_id = os.path.basename(fasta).split('.')[0]  # CJJ get prefix e.g. '4471'
        seqs = SeqIO.parse(fasta, 'fasta')

        if list_of_internal_outgroups:
            alignment = AlignIO.read(fasta, 'fasta')
            # Create an MultipleSeqAlignment object for ingroup sequences only:
            alignment_ingroup_seqs_only = []
            # Take first 10 sequences from the ingroup alignment for distance matrix calculations:
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

        # Populate the internal_outgroup_dict (potentially more than one sequence per taxon):
        for seq in seqs:
            seq_name_prefix = seq.name.split('.')[0]  # CJJ this assumes that there are no other dots ('.') in the
            # sequence name
            all_paralog_taxon_names.add(seq_name_prefix)
            if list_of_internal_outgroups and seq_name_prefix in list_of_internal_outgroups:
                internal_outgroup_dict[gene_id][seq_name_prefix].append(seq)

        # Filter internal outgroups with paralogs to retain the single sequence most divergent from the ingroup
        # sample:
        if list_of_internal_outgroups:
            logger.info(f'Detecting paralogs in internal outgroup sequences for gene {gene_id}...')
            internal_outgroup_dict_filtered = filter_internal_outgroups(internal_outgroup_dict,
                                                                        alignment_ingroup_only,
                                                                        gene_id)
            # Reassign internal_outgroup_dict[gene_id] dictionary entry to selected sequences:
            internal_outgroup_dict[gene_id] = internal_outgroup_dict_filtered[gene_id]

    # Read in external outgroups file, and create a dictionary of gene_id:list_of_seq_names, either for all seqs if
    # no external outgroup taxa specified, or for specified taxa only:
    external_outgroup_dict = defaultdict(list)
    if file_of_external_outgroups:  # CJJ dict not created if no outgroups file provided
        seqs = SeqIO.parse(file_of_external_outgroups, 'fasta')
        for seq in seqs:
            gene_id = seq.name.split('-')[-1]  # CJJ get gene id e.g. '4471'
            taxon = '-'.join(seq.name.split('-')[:-1])  # CJJ e.g. 'AMBTR'
            if list_of_external_outgroups:
                if taxon in list_of_external_outgroups:
                    seq.name = taxon  # CJJ i.e. we don't want the suffix e.g. '4471' present in the fasta file
                    seq.id = taxon
                    external_outgroup_dict[gene_id].append(seq)
            else:
                seq.name = taxon  # CJJ i.e. we don't want the suffix e.g. '4471' present in the fasta file
                seq.id = taxon
                external_outgroup_dict[gene_id].append(seq)

    all_external_outgroup_taxon_names = set([seq.name for gene_id, seq_list in external_outgroup_dict.items() for seq
                                             in seq_list])

    # Read in QC-d paralog files, add outgroup seqs, and write new fasta files ready for alignment:
    for fasta in glob.glob(f'{folder_of_qc_paralog_files}/*.selected.fa'):
        gene_id = re.sub('_[1-9]$', '', str(os.path.basename(fasta).split('.')[0]))
        gene_id_with_subtree_number = os.path.basename(fasta).split('.')[0]  # CJJ TEST

        seqs = list(SeqIO.parse(fasta, 'fasta'))
        if list_of_internal_outgroups:
            seqs_to_write = [seq for seq in seqs if seq.name.split('.')[0] not in list_of_internal_outgroups]
        else:
            seqs_to_write = seqs
        external_outgroup_seqs = external_outgroup_dict[gene_id]
        for taxon, sequence_list in internal_outgroup_dict[gene_id].items():  # add internal outgroup seqs
            seqs_to_write.extend(sequence_list)

        # seqs.extend(external_outgroup_seqs)  # add external outgroup seqs
        seqs_to_write.extend(external_outgroup_seqs)  # add external outgroup seqs

        # Write new files with outgroup sequences added (in the same directory as QC-d paralog files): #TODO: change
        #  CJJ this or it breaks Nextflow resume!!! CJJ is this done 19July2021?
        with open(f'{output_folder}/{gene_id_with_subtree_number}.outgroup_added.fasta', 'w') as outgroup_added:
            SeqIO.write(seqs_to_write, outgroup_added, 'fasta')

    # Write the IN and OUT taxon text file required by some paralogy resolution methods (MO, RT):
    if list_of_internal_outgroups:
        ingroup_taxon_names = [name for name in all_paralog_taxon_names if name not in list_of_internal_outgroups]
        logger.info(f'ingroup_taxon_names: {ingroup_taxon_names}')
    else:
        ingroup_taxon_names = [name for name in all_paralog_taxon_names]
        logger.info(f'ingroup_taxon_names: {ingroup_taxon_names}')
    with open(f'in_and_outgroups_list_{os.path.basename(folder_of_qc_paralog_files)}.txt', 'w') as group_list:
        if list_of_internal_outgroups:
            for taxon in list_of_internal_outgroups:
                group_list.write(f'OUT\t{taxon}\n')
        for taxon in all_external_outgroup_taxon_names:
            group_list.write(f'OUT\t{taxon}\n')
        for taxon in ingroup_taxon_names:
            group_list.write(f'IN\t{taxon}\n')

    return output_folder


def main(args):
    """
    Entry point for the resolve_paralogs.py script

    :param args: argparse namespace with subparser options for function main()
    :return:
    """

    # Initialise logger:
    logger = utils.setup_logger(__name__, 'logs_resolve_paralogs/08_align_selected_and_tree')

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

    # # Add outgroup sequences, both internal (if removed by the tree QC steps) and external (if a fasta file of
    # # external outgroup sequences is provided).
    # outgroups_added_folder = add_outgroup_seqs(results.original_paralog_gene_fasta_directory,
    #                                            results.gene_fasta_directory,
    #                                            results.internal_outgroups,
    #                                            results.external_outgroups_file,
    #                                            list_of_external_outgroups=results.external_outgroups)
    #
    # if not results.no_supercontigs:  # i.e. if it's a standard run.
    #     alignments_output_folder = mafft_align_multiprocessing(outgroups_added_folder,
    #                                                            algorithm=results.mafft_algorithm,
    #                                                            pool_threads=results.threads_pool,
    #                                                            mafft_threads=results.threads_mafft,
    #                                                            no_supercontigs=results.no_supercontigs,
    #                                                            use_muscle=results.use_muscle)
    #
    #     # Generate trees:
    #     if results.use_fasttree:
    #         fasttree_multiprocessing(alignments_output_folder,
    #                                  pool_threads=results.threads_pool,
    #                                  bootstraps=results.generate_bootstraps)  # Uses OpenMP with max threads default
    #     else:
    #         iqtree_multiprocessing(alignments_output_folder,
    #                                pool_threads=results.threads_pool,
    #                                iqtree_threads=results.threads_mafft,
    #                                bootstraps=results.generate_bootstraps)
    #
    # elif results.no_supercontigs:  # re-align with Clustal Omega.
    #     alignments_output_folder = mafft_align_multiprocessing(outgroups_added_folder,
    #                                                            algorithm=results.mafft_algorithm,
    #                                                            pool_threads=results.threads_pool,
    #                                                            mafft_threads=results.threads_mafft,
    #                                                            no_supercontigs=results.no_supercontigs,
    #                                                            use_muscle=results.use_muscle)
    #
    #     clustal_alignment_output_folder = clustalo_align_multiprocessing(alignments_output_folder,
    #                                                                      pool_threads=results.threads_pool,
    #                                                                      clustalo_threads=results.threads_mafft)
    #     # Generate trees:
    #     if results.use_fasttree:
    #         fasttree_multiprocessing(clustal_alignment_output_folder,
    #                                  pool_threads=results.threads_pool,
    #                                  bootstraps=results.generate_bootstraps)  # Uses OpenMP with max threads default
    #     else:
    #         iqtree_multiprocessing(clustal_alignment_output_folder,
    #                                pool_threads=results.threads_pool,
    #                                iqtree_threads=results.threads_mafft,
    #                                bootstraps=results.generate_bootstraps)

