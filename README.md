# HoloSimulator

**HoloSimulator** is a tool for generating synthetic short-read sequencing reads of eukaryotic hosts and associated microorganisms based on reference genomes and transcriptomes. 

## Installation

The easiest way to install **Holosimulator** is to create a conda environment containing the software and all its dependencies.

```sh
conda env create -f https://raw.githubusercontent.com/alberdilab/holosimulator/refs/heads/main/environment.yaml
```

## Usage

**HoloSimulator** has three main usage modes.

### Genomics

Produces genomic reads from reference genome sequences. Simulated reads are randomly generated for each genome, yet following the absolute or relative proportion of each genome within the sample.

```{sh}
holosimulator genomics -i input.tsv
```

### Transcriptomics

Produces transcriptomic reads from reference transcriptome sequences. The main difference from genomics is that each contig is treated as a transcript, and an expression pattern following a log-normal distribution is simulated instead of randomly assigning simulated reads to each contig.

```{sh}
holosimulator transcriptomics -i input.tsv
```

### Mutations

This operating mode has a completely different purpose to the above two. It introduces random mutations into a reference genome to reach a desired level of dissimilarity from the original sequence. This tool is useful to simulate variability within host or microbial genomes before simulating the sequencing reads. 

```{sh}
holosimulator mutations -i reference_genome.fa.gz -o manipulated_genome.fa.gz
```

## Reference sources

**HoloSimulator** allows using two types of genome and transcriptome sources. **Local files** are directly sourced from the environment in which **HoloSimulator** is running. However, it is also possible to rely on **remote files**, which are automatically fetched from the internet. Remote files typically derive from ftp sources such as https://ftp.ncbi.nlm.nih.gov/genomes/.

## Operating modes

For **HoloSimulator** has three operating modes.

- **CSV input:** the user introduces a table that defines the number or proportion of reads desired for the different genomes

- **Reference input:**

- **Prebuilt:**


## Examples

