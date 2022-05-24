# Benchmark Analysis for Determining Isoforms by Nanopore Data
Here I present the programs and the steps for carrying out a small benchmark pre-analysis for determining full-length mRNA isoforms by Nanopore RNA-seq data.

---

## Table of Contents
- [Input data](#input-data)
  - [Read Length Distribution](#read-length-distribution)
- [Programs](#programs)
  - [Reference-guided algorithms](#reference-guided-algorithms)
    - [FLAIR](#flair)
    - [bambu](#bambu)
    - [LIQA](#liqa)
  - [Short-read algorithm](#short-read-algorithm)
    - [cufflinks](#cufflinks)
- [Preliminary Benchmark Analysis](#preliminary-benchmark-analysis)

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

In addition, I have performed an alignment of Illumina short reads and Nanopore long reads against the isoform reference in order to observe the differences in coverage and read length between the two methods.

For this purpose I used Minimap2 to map both reads, converted the resulting SAM files to BAM and sorted and indexed them using samtools.

![IGV alignments](https://user-images.githubusercontent.com/82102364/169590965-9e8f0be7-0268-4b83-9365-f547c1f22ddf.png)

In this way we can observe that the Illumina reads are much shorter, since a single read does not cover the entire area of the sequin, while in Nanopore it does, although the coverage of the Illumina sample is much higher than that of Nanopore (32902 vs. 706).


---

## Programs
As for the programmes used to carry out this small benchmark pre-analysis, I have chosen 3 of the most recent algorithms that allow the analysis of RNA-seq data using the Oxford Nanopore technology, as the *de novo* reconstruction algorithms are much more complex to set up and run than the reference-guided ones.
These can be divided into:
- Reference-guided algorithms:
  - [FLAIR](#flair)
  - [bambu](#bambu)
  - [LIQA](#liqa)

In addition to these 3 that allow the use of dRNA-seq data, I will also use [Cufflinks](#cufflinks) to highlight possible differences with isoform detection using programmes that use Illumina data as input. 


> **Note:** 
> All programs that featured a Docker/Singularity image or Conda package have been installed in this way to facilitate the replicability of the results obtained in this analysis.

---

### Reference-guided algorithms

#### FLAIR

[**FLAIR**](https://github.com/BrooksLabUCSC/flair) (Full-Length Alternative Isoform analysis of RNA) is an algorithm for the correction, isoform definition and alternative splicing analysis of noisy reads that has been used mainly for native RNA developed by [Tang et al. (2018)](https://www.biorxiv.org/content/early/2018/09/06/410183).
FLAIR uses multiple alignment steps and splice site filters to increase confidence in the set of isoforms defined from noisy data. FLAIR was designed to be able to sense subtle splicing changes in nanopore data.

![flair](https://user-images.githubusercontent.com/82102364/169703545-22651da2-3a6d-40f9-8cea-3e9e6af711c2.png)


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

![bambu](https://user-images.githubusercontent.com/82102364/169703629-e21c910d-2b06-4fd7-93fb-378595e5d2b7.png)


To carry out the analysis using bambu, just run the [Rscript](https://github.com/3w1714n0/Benchmark_Analysis_full-length_isoforms/blob/main/RScripts/Bambu_Analysis.R).

> See the [bambu repository](https://github.com/GoekeLab/bambu) for further information.

Once the script is executed, three files are obtained:
- **counts_gene.txt**: a tab-delimited file with the counts per gene.
- **counts_transcript.txt**: a tab-delimited file with counts per transcript.
- **extended_annotations.gtf**: a tab-delimited file with the annotations made by the programme.

---

#### LIQA

[**LIQA**](https://github.com/WGLab/LIQA)(Long-read Isoform Quantification and Analysis) is an Expectation-Maximization based statistical method to quantify isoform expression and detect differential alternative splicing (DAS) events using long-read RNA-seq data. LIQA incorporates base-pair quality score and isoform-specific read length information to assign different weights across reads instead of summarizing isoform-specific read counts directly. Moreover, LIQA can detect DAS events between conditions using isoform usage estimates.

![Imagen1](https://user-images.githubusercontent.com/82102364/170011191-16c4c45b-7d18-4028-ae78-15f5160b0d4e.png)


To perform the analysis using LIQA, two steps must be carried out:

##### refgene

The first step is the transformation of isoforms into a compatible matrix based on a reference annotation file. In our case we will use the command for gtf files, although there is also an option for files in ucsc format.

```bash
liqa -task refgene -ref example.gtf -format gtf -out example.refgene
```
In our case:
- -ref: **rnasequin_annotation_2.4.gtf**

##### quantify

The next step is to quantify the expression of the isoforms. Although it is recommended to filter the bam file reads, this has not been done, as it has not been done with any of the above algorithms.

It is then necessary to sort and index the bam file, and give refgene_file, bam_file to LIQA to estimate isoform expression:

```bash
liqa -task quantify -refgene <refgene_file> -bam <bam_file> -out <output_file> -max_distance <max distance> -f_weight <weight of F function>
```
In our case:
- -refgene: **refgene_file** generated in the previous step.
- -bam: **flair.aligned.bam**, as the reads are aligned via minimap2 and we use the same input for all programmes.
- -max_distance: **10**
- -f_weight: **1**, as the example shown.

Finally, a file with the quantification of each isoform is obtained.

> See the [LIQA repository](https://github.com/WGLab/LIQA) for further information.

---

### Short-read algorithm

#### cufflinks
[**cufflinks**](https://github.com/cole-trapnell-lab/cufflinks) assembles transcripts, estimates their abundance and tests differential expression and regulation in RNA-Seq samples. It accepts aligned RNA-Seq reads and assembles the alignments into a parsimonious set of transcripts. Cufflinks then estimates the relative abundance of these transcripts based on the number of reads supported by each, taking into account biases in library preparation protocols.

As cufflinks only accepts aligned RNA-Seq reads, we must first align the sequence obtained during sequencing with the reference genome with the [minimap2 program](https://github.com/lh3/minimap2), using the following command:

```bash
minimap2 -a reference.fa reads.fq > alignment.sam
```

In our case: 
- reference.fa: **rnasequin_decoychr_2.4.fa**
- reads.fq: **k562+sequins_dRNA_albacore-2.1.3.fastq**

The next step is to change the format of the `aligment.sam` file into `aligment.bam` format, by means of the [samtools program](https://github.com/samtools/samtools), using the following command:

```bash
samtools view -S -b aligment.sam > alignment.bam
```

Next, we have to sort the `alignment.bam` file in genomic order with [samtools](https://github.com/samtools/samtools):

```bash
samtools sort alignment.bam -o alignment.sorted.bam
```

And finally, index the `alignment.sorted.bam` file:

```bash
samtools index alignment.sorted.bam
```

To facilitate the handling of Cufflinks and to apply a tool used during the degree, the following two steps have been carried out through [Galaxy Server](https://usegalaxy.org/), having available the workflow used.

---

## Preliminary Benchmark Analysis

To carry out the preliminary analysis of the reference analysis once we have managed to quantify at the isoform level with the three proposed programmes, I have generated an RScript that takes as input data from an auto-generated Excel file with the quantisation data for each programme and outputs a graph for each algorithm comparing the quantisation of that programme with the expected quantisation, available in the file `rnasequin_isoforms_2.4.tsv`.
