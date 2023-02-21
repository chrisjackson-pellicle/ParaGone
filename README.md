# ParaGone

Documentation in progress...

Current version: 0.0.3rc (Februrary 2023)

-----

### Purpose

ParaGone is a Python package (with external dependencies) designed to identify ortholog groups from a set of paralog sequences from multiple taxa. For phylogetic analyses, the correct inference of ortholog groups is critical in groups showing frequent gene or genome duplication, such as many families of land plants in which polyploidy is prevalent. Phylogenetic analysis of unrecognized paralogous gene copies can produce incorrect topologies, as the evolutionary history of gene families interferes with the evolutionary history of species lineages. 

ParaGone employs algorithmic approaches using gene tree topologies to infer ortholog groups. It makes use of the **paralogy resolution** (also described as **orthology inference**) approaches described and implemented by Yang and Smith 2014 [here][1]. These approaches have since been adapted for target capture datasets as described in the manuscript [here][2]. The original documentation and scripts can be found [here][3]. ParaGone provides a single installable package that uses many of the original Yang and Smith scripts, together with new pipeline steps and modules. The Yang and Smith have been updated and extended to include many new report and logging files, allowing the user to trace the processing of input paralogs fasta sequences through thge pipeline for all orthology inference algorithms. 


---

# Dependencies
* [Python](https://www.python.org/downloads/) 3.6 or later, along with the Python libraries:
    * [progressbar2](https://github.com/WoLpH/python-progressbar). The conda install can be found [here](https://anaconda.org/conda-forge/progressbar2).
    * [ete3](http://etetoolkit.org/)
    * [biopython](http://biopython.org/wiki/Main_Page) 1.79 or later
    
* [MAFFT](link)
* [ClustalOmega](link) 
* [HmmCleaner.pl](https://metacpan.org/dist/Bio-MUST-Apps-HmmCleaner/view/bin/HmmCleaner.pl)
* [Trimal](link)
* [IQTREE](link)  
* [FastTree](link). The conda install can be found [here](link).

---
# Setup

We strongly recommend installing ParaGone using [conda](https://docs.conda.io/en/latest/miniconda.html) with a new environment. This will install ParaGone, all required Python packages, and all required external programs. If you have conda installed and the channels `bioconda` and `conda-forge` have already been added, this can be done using the command:

```
conda create -n paragone -c chrisjackson-pellicle paragone
```

...followed by:

```
conda activate paragone
```

For full installation instructions, please see the wiki page:

https://github.com/chrisjackson-pellicle/ParaGone/wiki/Installation


Once all dependencies are installed you can download a test dataset [here](link), extract the `*.tar.gz` file, and execute the `run_paragone_test_dataset.sh` script from the `test_dataset` directory for a demonstration of ParaGone.


----

# Pipeline Input

Full instructions on running the pipeline, including a step-by-step tutorial using a small test dataset, are available on our wiki:

https://github.com/chrisjackson-pellicle/ParaGone/wiki/Tutorial

### Paralog sequences

A folder containing a fasta file for each gene. Each fasta file contains the 'main' contig selected by HybPiper for each sample. Where HybPiper has detected putative paralog contigs, these sequences are also included; in such cases, the main contig has the fasta header suffix `.main`, whereas putative paralogs have the suffix `.0`, `.1` etc.

### Outgroup sequences

Some of the paralogy resolution methods used in this pipeline require an outgroup sequence for each of your genes. These outgroup sequences can be provided in two ways.

1) Designating one or more taxa in your HybPiper paralog files as outgroups, via the `--internal_outgroup` parameter. For example, if your paralog `fasta` files contain sequences from the taxa `79686` and `79689`, you could designate these sequences as outgroups using `--internal_outgroup 79686 --internal_outgroup 79689`. 


2) Providing a fasta file (e.g. `outgroups.fasta`) containing 'external' outgroup sequences via the option `--external_outgroups_file outgroups.fasta`. The sequences in the file should have the same fasta header formatting and gene names as your HybPiper target file. For example, if you have used the Angiosperms353 target file for you HybPiper analysis, and you wish to use sequences from *Sesame* as your outgroup, your `outgroups.fasta` file might contain the following:

       >sesame-6995
       gtgggatatgaacaaaatccattgagcttgtattactgtta...
       >sesame-4757
       ctggtgcgtcgagcacttctcttgagcagcaacaatggcgg...
       >sesame-6933
       gaagtagatgctgtggtggtggaagcattcgacatatgcac...
    
       ...etc
    
Again, note that the gene identifier following the dash in the fasta headers (e.g. `6995` for header `>sesame-6995`) needs to correspond to a gene identifier in your target file. 

It's fine if your `outgroups.fasta` file contains additional sequences that you don't want to use. When running the pipeline (see below) you can optionally provide one or more taxon names using the parameter `--external_outgroup`, e.g. `--external_outgroup <taxon_id1> --external_outgroup <taxon_id2>`, and only these taxa will be included as outgroups. If this option isn't provided, all taxa/sequences in the `outgroups.fasta` file will be used.

**NOTE:** at a minimum, you must provide either 'internal' ingroups via the `--internal_outgroup` parameter, OR a file of 'external' outgroup sequences via the `--external_outgroups_file` parameter. You can also provide both 'internal' and 'external' outgroups.

----

# Pipeline Output

**IMPORTANT:** Unless you use the flag `--keep_intermediate_files` when running the command `paragone full_pipeline` (if running the pipeline with a single command) or `paragone final_alignments` (if running the pipeline in stages), most intermediate files and folders generated by the pipeline will be deleted. The information below assumes you have used the `--keep_intermediate_files` flag.

Running the pipeline will produce up to 28 output folders. The data in these folders can be useful for tracking and debugging. However, if you're just after the aligned `.fasta` sequences for each of your target genes (as output by each of the paralogy resolution methods), the main output folders of interest are probably:

    23_MO_final_alignments
    24_MI_final_alignments
    25_RT_final_alignments
    26_MO_final_alignments_trimmed
    27_MI_final_alignments_trimmed
    28_RT_final_alignments_trimmed

For a full description of ParaGone output, [see the wiki](https://github.com/chrisjackson-pellicle/ParaGone/wiki/Results-and-output-files).


-----
# Changelog

**0.0.1rc** *January, 2023*



[1]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4209138/ "Link to the Yang and Smith 2014 manuscript"
[2]: https://pubmed.ncbi.nlm.nih.gov/33978764/ "Link to Morales-Briones 2021 manuscript"
[3]: https://bitbucket.org/dfmoralesb/target_enrichment_orthology/src/master/ "Link to Yang and Smith Bitbucket"
