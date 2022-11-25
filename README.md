# ParaGone

Documentation in progress...

Current version: 1.0.1 (November 2022)

### Purpose

ToDo


#### Original Yang and Smith orthology inference manuscript, documentation and scripts


ParalogGrouper makes use of the **paralogy resolution** (also described as **orthology inference**) approaches described and implemented by Yang and Smith 2014 [here][1]. These approaches have since been adapted for target capture datasets as described in the bioRxiv manuscript [here][2]. The original documentation and scripts can be found [here][3].

---

# Dependencies
* [Python](https://www.python.org/downloads/) 3.6 or later, along with the Python libraries:
    * [progressbar2](https://github.com/WoLpH/python-progressbar). The conda install can be found [here](https://anaconda.org/conda-forge/progressbar2).
    * [ete3](http://etetoolkit.org/)
    * [biopython](http://biopython.org/wiki/Main_Page) 1.79 or later
    
* [MAFFT](link)
* [MUSCLE](link)
* [ClustalOmega](link) 
* [HmmCleaner.pl](https://metacpan.org/dist/Bio-MUST-Apps-HmmCleaner/view/bin/HmmCleaner.pl)
* [Trimal](link)
* [IQTREE](link)  
* [FastTree](link). The conda install can be found [here](link).

---
# Setup

We strongly recommend installing ParalogGrouper using [conda](https://docs.conda.io/en/latest/miniconda.html) with a new environment. This will install ParalogGrouper, all required Python packages, and all required external programs. If you have conda installed and the channels `bioconda` and `conda-forge` have already been added, this can be done using the command:

```
conda create -n paralog_grouper -c chrisjackson-pellicle paralog_grouper
```

...followed by:

```
conda activate paralog_grouper
```

For full installation instructions, please see our wiki page:

[address](link)


Once all dependencies are installed you can download a test dataset [here](link), extract the `*.tar.gz` file, and execute the `run_paralog_grouper_test_dataset.sh` script from the `test_dataset` directory for a demonstration of ParalogGrouper.


----

# Pipeline Input

Full instructions on running the pipeline, including a step-by-step tutorial using a small test dataset, are available on our wiki:

[address](link)

### Paralog sequences

A folder containing a fasta file for each gene. Each fasta file contains the 'main' contig selected by HybPiper for each sample. Where HybPiper has detected putative paralog contigs, these sequences are also included; in such cases, the main contig has the fasta header suffix .main, whereas putative paralogs have the suffix .0, .1 etc.

### Outgroup sequences

Some of the paralogy resolution methods used in this pipeline require an outgroup sequence for each of your genes. These outgroup sequences can be provided in two ways.

1) Designating one or more taxa in your HybPiper paralog files as outgroups, via the `--internal_outgroups <taxon1,taxon2,taxon3...>` option. For example, if your paralog `fasta` files contain sequences from the taxa `79686` and `79689`, you could designate these sequences as outgroups via `--internal_outgroups 79686,79689`. 


2) Providing a fasta file (e.g. `outgroups.fasta`) containing 'external' outgroup sequences via the option `--external_outgroups_file outgroups.fasta`. The sequences in the file should have the same fasta header formatting and gene names as your HybPiper target file. For example, if you have used the Angiosperms353 target file for you HybPiper analysis, and you wish to use sequences from *Sesame* as your outgroup, your `outgroups.fasta` file might contain the following:

       >sesame-6995
       gtgggatatgaacaaaatccattgagcttgtattactgtta...
       >sesame-4757
       ctggtgcgtcgagcacttctcttgagcagcaacaatggcgg...
       >sesame-6933
       gaagtagatgctgtggtggtggaagcattcgacatatgcac...
    
       ...etc
    
Again, note that the gene identifier following the dash in the fasta headers (e.g. '6995' for header '>sesame-6995') needs to correspond to a gene identifier in your target file. 

It's fine if your `outgroups.fasta` file contains additional sequences. When running the pipeline (see below) you can optionally provide one or more taxon names using the parameter `--external_outgroups <taxon1,taxon2,taxon3...>`, e.g. `--outgroups sesame`, and only these taxa will be included as outgroups. If this option isn't provided, all taxa/sequences in the `outgroups.fasta` file will be used. You can provide more than one outgroup taxon name using a comma-separated list, e.g. `--external_outgroups sesame,taxon2,taxon3` etc.

**NOTE:** at a minimum, you must provide either 'internal' ingroups via the `--internal_outgroups <taxon1,taxon2,taxon3...>` option, or a file of 'external' outgroup sequences via the `--external_outgroups_file outgroups.fasta` option.

----

# Pipeline Output

After running the pipeline, output can be found in the folder `results` (unless you have changed the name of the default output folder using the `--outdir <name>` parameter. This will consist of 23 subfolders. If you're just after the aligned .fasta files for each of your target genes as output by each of the paralogy resolution methods, the three main output folders of interest are probably:

- xx
- xx
- xx


For a full description of ParalogGrouper output, [see the wiki](link).


-----
# Changelog

**1.0.0 release candidate** *November, 2022*



[1]: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4209138/ "Link to the Yang and Smith 2014 manuscript"
[2]: https://www.biorxiv.org/content/10.1101/2020.08.21.261925v2 "Link to Yang 2021 bioarchives manuscript"
[3]: https://bitbucket.org/dfmoralesb/target_enrichment_orthology/src/master/ "Link to Yang and Smith Bitbucket"
