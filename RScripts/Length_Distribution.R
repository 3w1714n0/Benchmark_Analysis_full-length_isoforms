### Length Distribution
### Emiliano Navarro, 2022

# Clean workspace
rm(list= ls())

# Set working directory
setwd("C:/Users/enava/OneDrive - Universitat Autònoma de Barcelona/CUARTO/TFG/Análisis/sequins")

#1. Read data
reads_Ilum<-read.csv(file="ilum_length_distribution.txt", sep="", header=FALSE)
reads_Nano<-read.csv(file="nano_length_distribution.txt", sep="", header=FALSE)

#Nanopore mean length reads
weighted.mean(reads_Nano$V1, reads_Nano$V2)

#Illumina total reads
reads_Ilum$V2

#Nanopore total reads
sum(reads_Nano$V2)

# Adding sequentiation method
reads_Ilum$Type<-"Illumina"
reads_Nano$Type<-"Nanopore"

#Combinig the datasets
reads<-rbind(reads_Ilum, reads_Nano)

# 2. Plotting the length distribution
library(ggplot2)

ggplot(reads, mapping = aes(x = reorder(V1, V1), y = V2, fill = Type))+
  geom_bar(stat = "identity")+
  theme_bw()+
  scale_fill_manual(values = c("#fcc100","#2cc0e5"))+
  labs(x = "Read length", y="Number of reads")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  facet_grid(row = vars(Type), scales = "free")

# 1. Read data
ld <- read.csv(file="length_distribution.txt", sep="", header=FALSE)

# 2. Plotting the length distribution
library(ggplot2)

ggplot(ld, mapping = aes(x = order(as.integer(V2)), y = V3, fill = V4))+
  geom_bar(stat = "identity", width = 4)+
  theme_bw()+
  scale_fill_manual(values = c("#d40048","#2cc0e5"))+
  labs(x = "Read length", y="Number of reads")+
  facet_grid(row = vars(V4), scales = "free")+
  theme(legend.position = "none")+
  xlim(0,2200)
