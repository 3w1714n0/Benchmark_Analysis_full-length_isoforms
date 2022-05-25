### Benchmarking analysis
### Emiliano Navarro, 2022

# Clean workspace
rm(list= ls())

# Set working directory
setwd("C:/Users/enava/OneDrive - Universitat Autònoma de Barcelona/CUARTO/TFG/Análisis")

# Load the require library
library(readxl)

# Read the data
info <- read_excel("C:/Users/enava/OneDrive - Universitat Autònoma de Barcelona/CUARTO/TFG/Análisis/Counts per transcript.xlsx")

# Scatter Plots
library("ggpubr")

pdf("Scatterplots_benchmarking.pdf", height=6, width=6)

# Raw data Flair
ggscatter(info, x = "MIX_A", y = "Flair", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mix_A_isoform_counts", ylab = "Flair_isoform_counts",
          color = "#A80000")

# Raw data Bambu
ggscatter(info, x = "MIX_A", y = "Bambu", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mix_A_isoform_counts", ylab = "Bambu_isoform_counts",
          color = "#5A774E")

# Raw data LIQA
ggscatter(info, x = "MIX_A", y = "LIQA", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mix_A_isoform_counts", ylab = "LIQA_isoform_counts",
          color = "#002496")

# Raw data HTSeq
ggscatter(info, x = "MIX_A", y = "HTSeq", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Mix_A_isoform_counts", ylab = "HTSeq_isoform_counts",
          color = "#EEB500")
dev.off()
