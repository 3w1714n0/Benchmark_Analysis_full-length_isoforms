### Scimago Ranking Plot
### Emiliano Navarro, 2022

# Clean workspace
rm(list= ls())

# Set working directory
setwd("C:/Users/enava/Desktop")

# Read the data
journals <- c("Genome Biology", "Nucleic Acids Research", "Genome Research", "Genomics, Proteomics and Bioinformatics", "BMC Bioinformatics")
sjr <- c(9.371, 8.241, 5.707, 2.246, 1.246)
info <- data.frame(journals, sjr)

# Load the required library
library(ggplot2)

# Plot
ggplot(data = info, mapping = aes(x = reorder(journals, -sjr), y = sjr), fill = journals)+
  geom_col(fill = "#5A774E")+
  theme_bw()+
  labs(x = "Journals", y = "SJR (SCImago Journal Rank indicator)")
