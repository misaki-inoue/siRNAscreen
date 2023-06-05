library(dplyr)
setwd("C:/Users/misak/Desktop/ICR/siRNAscreen/siRNA_screen_validation_scripts/datafiles")

all_z <- read.csv("all_zscores.txt", sep = "\t", header = TRUE)

colnames(all_z) <- c("gene", "plate", "well", "zscore1", "zscore2", "zscore3")
myzscores<- all_z %>% filter(!is.na(gene))

library(stringr)
topwells <- myzscores %>% filter(str_detect(well, "^B|^C|^D"))
  #select rows with well name that starts with either B, C or D
botwells <- myzscores %>% filter(str_detect(well, "^I|^J|^K"))
colnames(botwells) <- c("gene", "plate", "well", "zscore4", "zscore5", "zscore6")

# Separate plate 1 and 2 before joining
topz_mut <- topwells %>% filter(plate == 1)
topz_res <- topwells %>% filter(plate == 2)
botz_mut <- botwells %>% filter(plate == 1) %>% select(gene, zscore4, zscore5, zscore6)
botz_res <- botwells %>% filter(plate == 2) %>% select(gene, zscore4, zscore5, zscore6)

# Assign different names for each individual siRNA
topz_mut <- topz_mut %>% group_by(gene) %>% 
  mutate(gene_id = paste(gene, row_number(), sep = "_"))
topz_mut <- topz_mut %>% ungroup() %>% select(gene_id, zscore1, zscore2, zscore3) %>%
  rename(gene = "gene_id")

topz_res <- topz_res %>% group_by(gene) %>% 
  mutate(gene_id = paste(gene, row_number(), sep = "_"))
topz_res <- topz_res %>% ungroup() %>% select(gene_id, zscore1, zscore2, zscore3) %>%
  rename(gene = "gene_id")

botz_mut <- botz_mut %>% group_by(gene) %>% 
  mutate(gene_id = paste(gene, row_number(), sep = "_"))
botz_mut <- botz_mut %>% ungroup() %>% select(gene_id, zscore4, zscore5, zscore6) %>%
  rename(gene = "gene_id")

botz_res <- botz_res %>% group_by(gene) %>% 
  mutate(gene_id = paste(gene, row_number(), sep = "_"))
botz_res <- botz_res %>% ungroup() %>% select(gene_id, zscore4, zscore5, zscore6) %>%
  rename(gene = "gene_id")


# Join top well zscores and bottom well zscores to make dataframe with 6 zscores per gene
six_z_mut <- left_join(topz_mut, botz_mut, by = "gene")
six_z_res <- left_join(topz_res, botz_res, by = "gene")

# Combine mutant and restored 6 z-scores
six_z_mut <- six_z_mut %>% mutate(status = "Mutant")
six_z_res <- six_z_res %>% mutate(status = "Restored")
six_zscores <- bind_rows(six_z_mut, six_z_res)

# Calculate average and standard deviation
six_zscores <- six_zscores %>% rowwise() %>%
  mutate(avg_z = mean(c(zscore1,zscore2,zscore3,zscore4,zscore5,zscore6))) %>%
  mutate(stdv = sd(c(zscore1,zscore2,zscore3,zscore4,zscore5,zscore6)))


# Export data
write.table(six_zscores, file = "six_zscores.txt", sep = "\t")


