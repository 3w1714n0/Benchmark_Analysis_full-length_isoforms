# Benchmark analysis for Nano-determining full-length isoforms
Here I present the programs and the steps for carrying out a benchmark analysis for determining full-length mRNA isoforms by Nanopore RNA-seq data.

---

## Table of Contents
- [Input data](#input-data)
  - [Read Length Distribution](#read-length-distribution)
- [Programs](#programs)
  - [flair](#flair)
  - [bambu](#bambu)
  - [CARNAC-LR](#carnac-lr)
  - [RATTLE](#rattle)
  - [cufflinks](#cufflinks)

---

## Input data

As for the datasets to perform a pre-screening study for this benchmark analysis, I have decided to use ***sequins***, synthetic RNA standards, available thanks to the [Garvan Institute of Medical Research](https://www.sequinstandards.com/).

For this initial analysis, I downloaded the following datasets:
- **Reference files**:
  - [rnasequin_annotation_2.4.gtf](https://sequins.s3.amazonaws.com/website/rna_resources/rnasequin_annotation_2.4.gtf) – annotation of sequin genes/isoforms on the decoy chromosome (120kb).
  - [rnasequin_sequences_2.4.fa](https://sequins.s3.amazonaws.com/website/rna_resources/rnasequin_sequences_2.4.fa) – sequences of all sequin isoforms (213kb).
- **Example libraries**:
  - [K562_SequinMixA.Rep2.R1.fq.gz (2.5Gb)](https://s3.amazonaws.com/sequins/website/libraries/K562_SequinMixA.Rep2.R1.fq.gz): Sequins (RNA  mixture A) was added to RNA extract from the K562 cell type. The sample then underwent library preparation (using the KAPA Stranded mRNA-Seq) and sequencing (using the **IIlumina HiSeq 2500**).
  - [k562sequins_dRNA_albacore2.1.3.tar.gz (33Mb)](https://s3.amazonaws.com/sequins/website/rna_resources/k562sequins_dRNA_albacore2.1.3.tar.gz): RNA sequins (Mixture A) was added to 500ng total RNA extract from the K562 cell type. The sample then underwent sequencing using **Nanopore MinION flowcell** (R9.4 chemistry), displayed 1323 active channels at commencement of run. Used sequencing protocol RNA001. Base called with Albacore version 2.1.3.

### Read Length Distribution

To highlight the differences between Illumina and Nanopore reads, I have plotted the read length distributions for each sequencing method, using the following command to obtain the read lengths:

```bash
awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} {for (l in lengths) {print l, lengths[l]}}' file.fastq >> file_with_read_lengths.txt
```
---

## Programs

### flair

### bambu

### CARNAC-LR

### RATTLE

### cufflinks
