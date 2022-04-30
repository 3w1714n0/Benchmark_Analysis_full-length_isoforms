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
  - [Short-read algorithm](#short-read-algorithm)
    - [cufflinks](#cufflinks)

---

## Input data

As for the datasets to perform a pre-screening study for this benchmark analysis, I have decided to use ***sequins***, synthetic RNA standards, available thanks to the [Garvan Institute of Medical Research](https://www.sequinstandards.com/).

For this initial analysis, I downloaded the following datasets:
- **Reference files**:
  - [rnasequin_decoychr_2.4.fa.gz](https://sequins.s3.amazonaws.com/website/rna_resources/rnasequin_decoychr_2.4.fa.gz) - decoy chromosome (chrIS) sequence the encodes all synthetic sequin gene loci (10Mb).
  - [rnasequin_annotation_2.4.gtf](https://sequins.s3.amazonaws.com/website/rna_resources/rnasequin_annotation_2.4.gtf) – annotation of sequin genes/isoforms on the decoy chromosome (120kb).
  - [rnasequin_sequences_2.4.fa](https://sequins.s3.amazonaws.com/website/rna_resources/rnasequin_sequences_2.4.fa) – sequences of all sequin isoforms (213kb).
  - [rnasequin_isoforms_2.4.tsv](https://sequins.s3.amazonaws.com/website/rna_resources/rnasequin_isoforms_2.4.tsv) – expected concentration (expression) of each sequin isoform in mixture (5kb).
- **Example libraries**:
  - [K562_SequinMixA.Rep2.R1.fq.gz](https://s3.amazonaws.com/sequins/website/libraries/K562_SequinMixA.Rep2.R1.fq.gz) (2.5Gb): Sequins (RNA  mixture A) was added to RNA extract from the K562 cell type. The sample then underwent library preparation (using the KAPA Stranded mRNA-Seq) and sequencing (using the **Illumina HiSeq 2500**).
  - [k562sequins_dRNA_albacore2.1.3.tar.gz](https://s3.amazonaws.com/sequins/website/rna_resources/k562sequins_dRNA_albacore2.1.3.tar.gz) (33Mb): RNA sequins (Mixture A) was added to 500ng total RNA extract from the K562 cell type. The sample then underwent sequencing using **Nanopore MinION flowcell** (R9.4 chemistry), displayed 1323 active channels at commencement of run. Used sequencing protocol RNA001. Base called with Albacore version 2.1.3.

### Read Length Distribution

To highlight the differences between Illumina and Oxford Nanopore reads, I have plotted the read length distributions for each sequencing method, using the following command to obtain the read lengths:

```bash
awk 'NR%4 == 2 {lengths[length($0)]++ ; counter++} END {for (l in lengths) {print l, lengths[l]}}' file.fastq >> file_with_read_lengths.txt
```

Once I obtained the length of the reads and their frequencies, using an [RScript](https://github.com/3w1714n0/Benchmark_Analysis_full-length_isoforms/blob/main/RScripts/Length_Distribution.R), I made the comparative graph to observe the differences, obtaining the following graph:

![Edited](https://user-images.githubusercontent.com/82102364/165273650-11073efd-c725-4a63-af37-ec31b0036f86.png)

We can observe that the difference in the average read lengths, since, while for Illumina we have a large number of reads of the same size (125 bp) for Nanopore we obtain a much lower performance, with fewer reads, but of a much greater length (917 bp), although it is true that the variability is very wide, obtaining reads ranging from 5 bp to ~3 kb, but for our preliminary analysis this is fine.

---

## Programs

As for the programmes used to carry out this benchmarking analysis, I have chosen 4 of the most recent algorithms that allow the analysis of RNA-seq data using Oxford Nanopore Technology.
These can be divided into:
- Reference-guided algorithms:
  - [flair](#flair)
  - [bambu](#bambu)
- *De novo* reconstruction algorithms:
  - [CARNAC-LR](#carnac-lr)   
  - [RATTLE](#rattle)

In addition to these 4 that allow the use of dRNA-seq data, I will also use [Cufflinks](#cufflinks) to highlight possible differences with isoform detection using programmes that use Illumina data as input. 


> **Note:** 
> All programs that featured a Docker/Singularity image or Conda package have been installed in this way to facilitate the replicability of the results obtained in this analysis.

---

### Reference-guided algorithms

#### flair

[FLAIR](https://github.com/BrooksLabUCSC/flair) (Full-Length Alternative Isoform analysis of RNA) is an algorithm for the correction, isoform definition and alternative splicing analysis of noisy reads that has been used mainly for native RNA developed by [Tang et al. (2018)](https://www.biorxiv.org/content/early/2018/09/06/410183).
FLAIR uses multiple alignment steps and splice site filters to increase confidence in the set of isoforms defined from noisy data. FLAIR was designed to be able to sense subtle splicing changes in nanopore data.
 <p align="center">
  <img src='https://people.ucsc.edu/~atang14/flair/flair_workflow_compartmentalized.png' alt='flair workflow' width='680'/>
</p>
 
As we can see in the image above, this algorithm has multiple modules to follow, so I will therefore show below the commands I have used to carry out the analysis, but you can visit the [FLAIR repository](https://github.com/BrooksLabUCSC/flair) for further information.

##### flair align

Align reads.

```bash
singularity exec flair_latest.sif python3 ./flair/flair.py align -g genome.fa -r <reads.fq.gz>|<reads.fq>|<reads.fa> -nvra -o output_aligned
```
In our case:
- -g: **rnasequin_decoychr_2.4.fa**
- -r: **k562+sequins_dRNA_albacore-2.1.3.fastq**

##### flair correct

Correct reads using an annotation file.

```bash
singularity exec flair_latest.sif python3 ./flair/flair.py correct -q output_aligned.bed -g genome.fa -f annotation.gtf -c chromsize.fa.fai -nvrna -o output_correct
```
In our case:
- -g: **rnasequin_decoychr_2.4.fa**
- -f: **rnasequin_annotation_2.4.gtf**
- -c: **rnasequin_decoychr_2.4.fa.fai**

##### flair collapse

Defines high-confidence isoforms from corrected reads.

```bash
singularity exec flair_latest.sif python3 ./flair/flair.py collapse -g genome.fa -r <reads.fq>|<reads.fa> -f annotation.gtf -q output_corrected.bed --trust_ends
```
In our case:
- -g: **rnasequin_decoychr_2.4.fa**
- -r: **k562+sequins_dRNA_albacore-2.1.3.fastq**
- -f: **rnasequin_annotation_2.4.gtf**

##### flair quantify

Convenience function to quantifying FLAIR isoform usage across samples using minimap2, creating a manifest.tsv tab-separated file with the following structure:

```tsv
sample1 conditionA  batch1  ./sample1_reads.fq
```
Command to quantifying:

```bash
singularity exec flair_latest.sif python3 ./flair/flair.py quantify -r reads_manifest.tsv -i flair.collapse.isoforms.fa --trust_ends
```

After all the analysis, we obtain a .gtf file with all the isoforms found by FLAIR and a table with their quantification, which will be later compared with the quantifications of the sequins and the predictions made by the other programmes.

---

#### bambu

[***bambu***](https://github.com/GoekeLab/bambu) is a R package for multi-sample transcript discovery and quantification using long read RNA-Seq data. You can use bambu after read alignment to obtain expression estimates for known and novel transcripts and genes. The output from bambu can directly be used for visualisation and downstream analysis such as differential gene expression or transcript usage.

To carry out the analysis using bambu, just run the [Rscript](https://github.com/3w1714n0/Benchmark_Analysis_full-length_isoforms/blob/main/RScripts/Bambu_Analysis.R).

> See the [bambu repository](https://github.com/GoekeLab/bambu) for further information.

Once the script is executed, three files are obtained
- **counts_gene.txt**: a tab-delimited file with the counts per gene.
- **counts_transcript.txt**: a tab-delimited file with counts per transcript.
- **extended_annotations.gtf**: a tab-delimited file with the annotations made by the programme.

---

### *De novo* reconstruction algorithms

#### CARNAC-LR

#### RATTLE

### Short-read algorithm

#### cufflinks
