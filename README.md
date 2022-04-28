# Benchmark Analysis for Determining Isoforms by Nanopore Data
Here I present the programs and the steps for carrying out a benchmark analysis for determining full-length mRNA isoforms by Nanopore RNA-seq data.

---

## Table of Contents
- [Input data](#input-data)
  - [Read Length Distribution](#read-length-distribution)
- [Programs](#programs)
  - [Reference-guided algorithms](#reference-guided-algorithms)
    - [flair](#flair)
    - [bambu](#bambu)
  - [*De novo* reconstruction algorithms](#de-novo-reconstruction-algorithms)
    - [CARNAC-LR](#carnac-lr)
    - [RATTLE](#rattle)
  - [Illumina algorithm](#illumina-algorithm)
    - [cufflinks](#cufflinks)

---

## Input data

As for the datasets to perform a pre-screening study for this benchmark analysis, I have decided to use ***sequins***, synthetic RNA standards, available thanks to the [Garvan Institute of Medical Research](https://www.sequinstandards.com/).

For this initial analysis, I downloaded the following datasets:
- **Reference files**:
  - [rnasequin_annotation_2.4.gtf](https://sequins.s3.amazonaws.com/website/rna_resources/rnasequin_annotation_2.4.gtf) – annotation of sequin genes/isoforms on the decoy chromosome (120kb).
  - [rnasequin_sequences_2.4.fa](https://sequins.s3.amazonaws.com/website/rna_resources/rnasequin_sequences_2.4.fa) – sequences of all sequin isoforms (213kb).
- **Example libraries**:
  - [K562_SequinMixA.Rep2.R1.fq.gz](https://s3.amazonaws.com/sequins/website/libraries/K562_SequinMixA.Rep2.R1.fq.gz) (2.5Gb): Sequins (RNA  mixture A) was added to RNA extract from the K562 cell type. The sample then underwent library preparation (using the KAPA Stranded mRNA-Seq) and sequencing (using the **Illumina HiSeq 2500**).
  - [k562sequins_dRNA_albacore2.1.3.tar.gz](https://s3.amazonaws.com/sequins/website/rna_resources/k562sequins_dRNA_albacore2.1.3.tar.gz) (33Mb): RNA sequins (Mixture A) was added to 500ng total RNA extract from the K562 cell type. The sample then underwent sequencing using **Nanopore MinION flowcell** (R9.4 chemistry), displayed 1323 active channels at commencement of run. Used sequencing protocol RNA001. Base called with Albacore version 2.1.3.

### Read Length Distribution

To highlight the differences between Illumina and Nanopore reads, I have plotted the read length distributions for each sequencing method, using the following command to obtain the read lengths:

```bash
awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths) {print l, lengths[l]}}' file.fastq >> file_with_read_lengths.txt
```

Once I obtained the length of the reads and their frequencies, using an [RScript](https://github.com/3w1714n0/Benchmark_Analysis_full-length_isoforms/blob/main/RScripts/Length_Distribution.R), I made the comparative graph to observe the differences, obtaining the following graph:

![Edited](https://user-images.githubusercontent.com/82102364/165273650-11073efd-c725-4a63-af37-ec31b0036f86.png)

We can observe that the difference in the average read lengths, since, while for Illumina we have a large number of reads of the same size (125 bp) for Nanopore we obtain a much lower performance, with fewer reads, but of a much greater length (917 bp), although it is true that the variability is very wide, obtaining reads ranging from 5 bp to ~3 kb, but for our preliminary analysis this is fine.

---

## Programs

As for the programmes used to carry out this benchmarking analysis, I have chosen 4 of the most recent algorithms that allow the analysis of RNA-seq data using Nanopore.
These can be divided into:
- Reference-guided algorithms:
  - [flair](#flair)
  - [bambu](#bambu)
- *De novo* reconstruction algorithms:
  - [CARNAC-LR](#carnac-lr)   
  - [RATTLE](#rattle)

In addition to these 4 that allow the use of dRNA-seq data, I will also use [Cufflinks](#cufflinks) to highlight possible differences with isoform detection using programmes that use Illumina data as input. 

### Reference-guided algorithms

#### flair

#### bambu

### *De novo* reconstruction algorithms

#### CARNAC-LR

#### RATTLE

### Illumina algorithm

#### cufflinks
