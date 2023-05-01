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
    