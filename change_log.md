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
    