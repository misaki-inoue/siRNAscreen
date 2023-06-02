library(dplyr)
setwd("C:/Users/misak/Desktop/ICR/siRNA/datafiles_4_12")

all_z <- read.csv("all_zscores.txt", sep = "\t", header = TRUE)

colnames(all_z) <- c("gene", "plate", "well", "zscore1", "zscore2", "zscore3")
myzscores<- all_z %>% filter(!is.na(gene))

library(stringr)
topwells <- myzscores %>% filter(str_detect(well, "^B|^C|^D|^E"))
  #select rows with well name that starts with either A, B, or c
botwells <- myzscores %>% filter(str_detect(well, "^I|^J|^K|^L"))
colnames(botwells) <- c("gene", "plate", "well", "zscore4", "zscore5", "zscore6")

# Separate plate 1 and 2 before joining
topz_mut <- topwells %>% filter(plate == 1)
topz_res <- topwells %>% filter(plate == 2)
botz_mut <- botwells %>% filter(plate == 1) %>% select(gene, zscore4, zscore5, zscore6)
botz_res <- botwells %>% filter(plate == 2) %>% select(gene, zscore4, zscore5, zscore6)

# Join top well zscores and bottom well zscores to make dataframe with 6 zscores per gene
six_z_mut <- left_join(topz_mut, botz_mut, by = "gene")
six_z_res <- left_join(topz_res, botz_res, by = "gene")

# Take average all six z-scores and extract that column
avg_z_mut <- six_z_mut %>% rowwise() %>% mutate(avg_z_mut = mean(
  c(zscore1,zscore2,zscore3,zscore4,zscore5,zscore6))) %>% select(gene, avg_z_mut)
avg_z_res <- six_z_res %>% rowwise() %>% mutate(avg_z_res = mean(
  c(zscore1,zscore2,zscore3,zscore4,zscore5,zscore6))) %>% select(gene, avg_z_res)

z_diff <- left_join(avg_z_mut, avg_z_res, by = "gene")
z_diff <- z_diff %>% rowwise() %>% mutate(mut_res_diff = avg_z_mut - avg_z_res)

# Export data
setwd("C:/Users/misak/Desktop/ICR/siRNA")
six_z_res_sel <- six_z_res %>% select(
  gene, zscore1, zscore2, zscore3, zscore4, zscore5, zscore6
)
all_six_z <- left_join(six_z_mut, six_z_res_sel, by = "gene", suffix = c("mut","res"))
  # 12 zscore columns, 6 mutant & 6 restored

write.table(all_six_z, file = "all_six_zscores.txt", sep = "\t")

write.table(z_diff, file = "zscores_difference.txt", sep = "\t")

