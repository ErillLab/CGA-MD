# CGA-MD
Comparative genomics-assisted motif discovery

This comparative genomics pipeline is based on the [CGA-MD pipeline](https://github.com/ErillLab/CGA-MD/) implemented by [Maria Navarro](https://github.com/manaves), which was in turn based on the [blast_seq_filter code](https://github.com/ErillLab/blast_seq_filter) codebase developed by [Issac Chaudry](https://github.com/ichaudr1).

## Basic idea
The basic idea in CGA-MD is to enhance motif discovery through comparative genomics. Whether it is from ChIP-seq or RNA-seq data, or simply from the suspicion that a particular gene (or set of genes) should be regulated by a given transcription factor (TF), a main problem with identifying the TF-binding motif is the lack of enough data to perform the inference (e.g. with MEME).

The comparative genomics approach makes the (not very strong) assumption that the TF-binding motif is conserved within a given taxonomic group, and that in said group there are enough orthologs of the proteins that we have evidence of being regulated in a given species. Based on this premise, the idea is very simple:

- identify orthologs within the taxonomic group for each of the proteins on which we have evidence of regulation by the TF
- identify genome records harboring the genes encoding these orthologs
- pull down the regions upstream of the all these ortholog-encoding genes
	- in doing so, pull down also the upstream regions of genes preceding the ortholog-encoding genes, in case they are in an operon (and we want to grab the operon lead gene, which is where the regulatory activity will be)
- perform motif discovery on all these sequences
- prioritize inferred motifs based on their presence in diverse sequences, motif structure...

# Configuration file guidelines

This document provides detailed information on the configuration parameters used in the project. The configuration file follows a schema defined by Cerberus, ensuring that all inputs meet the specified criteria. Below are the parameters and their validation rules.

## 1. Entrez Parameters

The `entrez_parameters` section configures parameters related to Entrez queries.

- **request_limit** (integer, required): the maximum number of retries to be attempted on NCBI server. Must be at least 1. Recommended: 5.
- **sleep_time** (number, required): the time (in seconds) to wait between requests. Must be a non-negative number. Recommended: 1.
- **email** (string, required): a valid email address used for Entrez queries. Must match the regex pattern for email addresses.
- **api_key** (string, optional): the API key for Entrez. Can be obtained with an NCBI user account. Can be `null`.

## 2. BLAST Parameters

The `blast_parameters` section configures parameters related to BLAST searches.

- **db** (string, optional): the database to search against (e.g. `nr`). Can be `null`. Recommended: `nr_cluster_seq`.
- **tax_IDs** (list, optional): a list of taxonomy IDs to constrain the BLAST search. Can be `null` (unconstrained).
- **e_value** (number, optional): the E-value threshold for the BLAST search. Can be `null`. Recommended: 10E-10.
- **query_coverage** (number, optional): the query coverage fraction for establishing orthology. Can be `null`. Recommended: 0.75.
- **max_hits** (integer, optional): the maximum number of BLAST hits to return. Must be at least 1 if specified. Can be `null`. Recommended: 500.

## 3. Sequences Parameters

The `sequences_parameters` section configures parameters related to sequence processing.

- **max_intergenic_size** (number, required): the maximum intergenic size. Intergenic regions larger than this will signal that operon lead gene has been reached. Cannot be `null`.
- **min_intergenic_size** (number, required): the minimum intergenic size. Intergenic regions smaller than this will be ignored (not retrieved). Cannot be `null`.
- **upstream_size_region** (number, required): the size of the upstream region to explore. This is the region upstream of the gene that will be analyzed for intergenic regions forming a potential operon. Cannot be `null`. Recommended: 10000.
- **downstream_size_region** (number, required): the size of the downstream region. This only ensures that part of the region downstream the gene end will be retrieved. Cannot be `null`. Recommended: 50.
- **max_sequence_length** (integer, optional): the maximum sequence length for retrieved intergenic regions. For regions over this limit, only the max_sequence_length upstream of the current gene will be retrieved. Can be `null`.
- **crawl_back** (boolean, required): determines whether the system will scan backwards from the gene, trying to reach a putative operon lead gene. Cannot be `null`.
- **fuse_intergenic** (boolean, required): determines whether intergenic sequences pertaining to a putative operon identified during crawl back will be concatenated or not for motif discovery. If `false`, each intergenic region is retrieved as an independent sequence for motif discovery. This impacts the filtering step ahead of motif discovery, considering, if `true`, the whole set of intergenic sequences in the operon as a single sequence. Cannot be `null`.


## 4. Maximum Identity

- **max_identity** (float, required): the maximum identity threshold. This will determine which intergenic sequences are discarded because they are too similar to other retrieved intergenic regions. Must be between 0.0 and 1.0 inclusive. Cannot be `null`. Recommended: 0.75.

## 5. MEME Parameters

The `meme_parameters` section configures parameters for MEME motif discovery.

- **mod** (string, required): the model to use. Can be ```oops|zoops|anr```. Cannot be `null`.
- **nmotifs** (integer, required): the number of motifs to find. Cannot be `null`. Recommended: 10.
- **minw** (integer, required): the minimum motif width. Must be at least 0. Cannot be `null`. Recommended: 10.
- **maxw** (integer, required): the maximum motif width. Cannot be `null`. Recommended: 22.
- **revcomp** (boolean, required): whether to search both strands. Cannot be `null`. Recommended: ```True```.
- **pal** (boolean, required): whether to search for palindromic motifs. Cannot be `null`.

## 6. Run Only MEME

- **run_only_meme** (boolean, required): flag to indicate if only MEME should be run. MEME then will be run on the existing `sequences_filtered.fasta` file. Used to rerun MEME and the downstream MEME analysis pipeline with updated paramters. Cannot be `null`.

## 7. Output Parameters

The `output_parameters` section configures parameters related to output directories.

- **folder_name** (string, required): the name of the output folder. Cannot be `null`.
- **folder_rerun_meme** (string, optional): the folder name for rerunning MEME. Can be `null`.

## 8. Input Records

- **input_records** (list, required): a list of input records (NCBI protein accession numbers). Must contain at least one record. Cannot be `null`.

# Execution Instructions

## Step-by-Step Guide

### 1. Open the Terminal

Ensure you have access to the command line interface on your operating system. This could be Terminal on macOS or Linux, or Command Prompt/PowerShell on Windows.

### 2. Navigate to the Project Directory

Change your current directory to the root directory of the project. This directory contains a subdirectory ```/scripts``` where the main script is located. You can use the `cd` command followed by the path to the project directory. For example:

```sh
cd path/to/your/project
```

### 3. Execute the Main Script

To run the main pipeline, use the Python interpreter followed by the script name (with ```scripts/``` prepended) and the required root directory argument. The root directory is the main directory where the pipeline will operate. The command structure is as follows:

```sh
python scripts/main.py <root_directory>
```

Your root directory must contain a config folder with the config.json file in it.

#### Example

If your root directory is located at `/run` within the folder project, the command would be:

```sh
python scripts/main.py run
```
