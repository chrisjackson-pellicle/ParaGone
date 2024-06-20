**1.1.0** *20 June, 2024*

- Add new dependency: Julia
- Replaced HmmCleaner with TAPER for sequence-specific cleaning/QC.
- Remove use of Bio.Align.Applications.MafftCommandline and Bio.Align.Applications.ClustalOmegaCommandline due to deprecation.
- Ensure an internal outgroup sequence is always selected for each taxon when multiple sequences are present (issue#5).
- Add check that internal outgroup sequences are found in the folder provided to command `align_selected_and_tree` (issue#6).
- Update error handling and logging for the steps of the pipeline that use multiprocessing.
- Update commands in `run_paragone_test_dataset.sh` to match output directory names.

**1.0.0** *25 August, 2023*

- Bugfix: change type from `int` to `float` for `--trimal_resoverlap` in all relevant subparsers.
- Add check that `--trimal_resoverlap` is provided if `--trimal_seqoverlap` is used, and vice-versa.
- Add module `version.py` for a single location of ParaGone version number.

**0.0.14** *7th July, 2023*

- Change the parameter `new_mo_algorithm` to `mo_algorithm_paragone` for subparsers `parser_full_pipeline` and `parser_prune_paralogs`

**0.0.13** *6th July, 2023*

- Remove the positional parameter `selected_alignment_directory` from subparser `add_align_selected_and_tree_parser`, as `align_selected_and_tree.py` now hardcodes this path as `09_sequences_from_qc_trees`

**0.0.12** *3rd July, 2023*

- Add full trimal options to steps that use trimal: `check_and_align`, `align_selected_and_tree`, and `final_alignments`, as well as `full_pipeline`. 
- Bugfix in `prune_paralogs_mo.py`: correct error when parsing outgroup tip labels in function `prune_paralogs_from_rerooted_homotree_cjj()`.

**0.0.11** *1st May, 2023*

- Following QC but before paralog pruning, trees are now parsed to produce two additional reports in *.tsv format: `per_locus_paralogy_report_post_tree_qc.tsv` and `per_taxon_paralogy_report_post_tree_qc.tsv`. 
- Change optional flag `--original_mo_algorithm` to `--new_mo_algorithm` for subcommand `prune_paralogs` and `full_pipeline`; if used, apply a re-implementation of the MO algorithm rather than the original Yang and Smith MO algorithm. Default is to apply the original Yang and Smith MO algorithm.

**0.0.10** *31st March, 2023*

- Fix bug in `prune_paralogs_mo.py`: update outgroup trimming regex in new MO implementation to allow for optional bootstrap support values, and/or tree lengths in scientific notation (e.g. `2.079e-06`)    


**0.0.9** *23rd March, 2023*

- Add the flag `--debug` to the `prune_paralogs` subcommand. When supplied, additional details will be written to the log file. This can increase the size of the log file substantially.
- Fix bug in `prune_paralogs_mo.py` where the default implementation of the MO algorithm incorrectly parsed trees rooted with a single outgroup taxon.

**0.0.8** *17th March, 2023*

- Fix bug in `prune_paralogs_rt.py` where ingroups with fewer than the minimum required taxa caused an endless loop.
- Add optional flag `--original_mo_algorithm` to subcommand `prune_paralogs` and `full_pipeline`; if used, apply the original Yang and Smith MO algorithm. Default is to apply a re-implementation of the MO algorithm.
- Prevent full MAFFT alignments from being written to log files on Linux.


**0.0.7** *9th March, 2023*

- Make the MAFFT option `--adjustdirection` (https://mafft.cbrc.jp/alignment/software/adjustdirection.html) optional for the first alignment step of the pipeline; default is not to use this flag.
- Fix bug in MO pruning step when outgroups not monophyletic.
- Rename flag `--no_stitched_contigs` to `--use_clustal`. Only run MAFFT first if flag `--mafft_adjustdirection` is also supplied.
- For the tree QC step involving removal of long tip-branches, replace the original Yang and Smith 2014 script `trim_tips.py` with a wrapper around `TreeShrink`.
- Remove the trim tips step from the MI pruning stage; seqs are extracted from the MI pruned trees, so long branches are not a concern at this point.
- HmmCleaner now runs in parallel using multiprocessing.
- Added subcommand `delete_intermediate_files`.
    