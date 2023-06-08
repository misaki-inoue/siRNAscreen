library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)

setwd("C:/Users/misak/Desktop/ICR/siRNAscreen/siRNA_screen_validation_scripts/datafiles")

six_zscores <- read.csv("six_zscores.txt", sep = "\t", header = TRUE)

# Separate and reshape dataframe
avg_sd <- six_zscores %>% select(gene, status, avg_z, stdv)
long_zscore <- six_zscores %>% select(gene, status, zscore1, zscore2, zscore3,
                                      zscore4, zscore5, zscore6)
long_zscore <- long_zscore %>% pivot_longer(cols = 3:8, names_to = "replicates",
                                            values_to = "zscore")

# Cleaning up long_zscore
long_zscore$status <- factor(long_zscore$status, levels = c("Restored", "Mutant"))
long_zscore <- long_zscore %>% select(!replicates)
long_zscore <- long_zscore %>% mutate(gene = if_else(
  gene == "DHODH_1", "DHODH", gene)) # change just DHODH_1 and keep others the same
long_zscore$gene <- factor(long_zscore$gene, levels = c(
  "CAD_1","CAD_2","CAD_3","CAD_4","UCP2_1","UCP2_2","UCP2_3","UCP2_4","DHODH"
))


# Calculate mean and standard deviation
long_zscore.summary <- long_zscore %>% group_by(gene, status) %>%
  summarize(sd = sd(zscore), zscore = mean(zscore)) 
  # sd must come before mean or else it will try to calculate sd of the mean which gives NA


# Plot data

p <- ggplot(long_zscore, aes(x=gene, y=zscore, color = status)) +
  geom_jitter(position=position_dodge(0.8)) +
  geom_crossbar(
    aes(ymin = zscore-sd, ymax = zscore+sd),
    data = long_zscore.summary, 
    width = 0.5,
    position = position_dodge(0.8)) +
  scale_color_manual(name = "BRCA1", values = c("#0080ff", "#ff6666")) 
p

p1 <- p + theme_classic() +
  ylab("Z-score") + xlab("Gene") + theme(text = element_text(size = 12)) +
  geom_hline(yintercept = 0, 
             linetype="dashed", 
             color = "grey10", 
             alpha=0.4,
             linewidth=0.5)
p1
