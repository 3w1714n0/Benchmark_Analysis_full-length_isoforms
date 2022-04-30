### Bambu Analysis
### Emiliano Navarro, 2022

# Clean workspace
rm(list= ls())

# Set working directory
setwd("C:/Users/enava/OneDrive - Universitat Autònoma de Barcelona/CUARTO/TFG/Análisis/Bambu")

# Load the package
library(bambu)

# Read the aligned reads (bam), the reference genome sequence (fa) and reference genome annotations (gtf) data
test.bam <- system.file("extdata", "flair.aligned.bam", package = "bambu")

fa.file <- system.file("extdata", "rnasequin_decoychr_2.4.fa", package = "bambu")

gtf.file <- system.file("extdata", "rnasequin_annotation_2.4.gtf", package = "bambu")

# Prepare the annotations
bambuAnnotations <- prepareAnnotations(gtf.file)

# Run the analysis and save it in a summarized experiment object
se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file)

# Write bambu outputs to files
writeBambuOutput(se, path = "C:/Users/enava/OneDrive - Universitat Autònoma de Barcelona/CUARTO/TFG/Análisis/Bambu")
