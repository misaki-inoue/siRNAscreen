setwd("C:/Users/misak/Desktop/ICR/siRNA")

library(dplyr)

zscores <- read.csv("./datafiles_4_12/zscores.txt", sep="\t", header=TRUE)
toptable <- read.csv("./datafiles_4_12/TopTable.txt", sep="\t", header=TRUE)
genelibrary <- read.csv("./datafiles_4_12/Genelibrary.txt", sep="\t", header=TRUE)


par_data <- toptable %>% filter(plate == 1 & wellAnno == "sample")

res_data <- toptable %>% filter(plate == 2)

par_data <- par_data %>% filter(wellAnno != "empty")
res_data <- res_data %>% filter(wellAnno != "empty")

# # make dataframe with well, gene, and 3 normalized medians 
# par_medians <- par_data %>% 
#   select(well, normalized_r1_ch1, normalized_r2_ch1, normalized_r3_ch1 )
#   # normalized_r1_ch1 is normalized intensity values... welp

par_zsc <- zscores %>% filter(plate == 1)
res_zsc <- zscores %>% filter(plate == 2)
  # res_zsc doesn't have any gene

gene_well <- par_zsc %>% select(gene, well)
res_zsc <- res_zsc %>% select(-gene)
res_zsc <- left_join(res_zsc, gene_well, by = "well")

par_zsc <- par_zsc %>% filter(!is.na(gene))
par_zsc <- par_zsc %>% distinct(gene, .keep_all = TRUE)
res_zsc <- res_zsc %>% filter(!is.na(gene))
res_zsc <- res_zsc %>% distinct(gene, .keep_all = TRUE)

two_zsc <- left_join(par_zsc, res_zsc, by = "gene", suffix = c(".p", ".r"))


library(ggplot2)

p <- ggplot(two_zsc, aes(x=zscore.p, y=zscore.r)) + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color="red", 
              linetype="dashed", linewidth=0.8) +
  geom_vline(xintercept = 0, linetype="solid", 
             color = "grey", linewidth=0.5) +
  geom_hline(yintercept=0, linetype="solid", color = "grey") +
  scale_x_continuous(breaks= seq(-10, 4, by = 2)) +
  scale_y_continuous(breaks = seq(-10, 4, by = 2)) +
  ylab("BRCA1 Restored z-scores") + xlab("BRCA1 Mutant z-scores")
p

two_zsc <- two_zsc %>% mutate(sum_z = abs(zscore.p) + abs(zscore.r))
top <- 15
top_genes <- bind_rows( 
  two_zsc %>% 
    filter(zscore.p < -1.96|zscore.r < -1.96) %>% 
    arrange(desc(abs(sum_z))) %>% 
    head(top)
)

library(ggrepel)
# adding label to data

p2 <- p +
  geom_label_repel(data = top_genes,
                   mapping = aes(zscore.p, zscore.r, label = gene),
                   size = 3.5, max.overlaps = 50, min.segment.length = 0)
p2

# histogram of luminescence distribution

library(tidyr)
samples <- toptable %>% filter(wellAnno == "sample") %>% 
  select(plate, raw_r1_ch1, raw_r2_ch1, raw_r3_ch1)

samples <- samples %>% gather(key = "replicate", value = "luminescence",
                              raw_r1_ch1:raw_r3_ch1)
samples <- samples %>% mutate(BRCA1 = case_when(plate == 1 ~ "Mutant",
                                                plate == 2 ~ "Restored"))
  # cannot use plate column for groups when making histogram because they are numeric
options(scipen = 999)
h <- ggplot(samples, aes(x=luminescence, fill=BRCA1, color=BRCA1)) +
  geom_histogram(position="identity", bins = 40, alpha=0.5) +
  ylab("Count") + xlab("Raw Luminescence")
h

samples_norm <- toptable %>% filter(wellAnno == "sample") %>% 
  select(plate, normalized_r1_ch1, normalized_r2_ch1, normalized_r3_ch1)
samples_norm <- samples_norm %>% gather(key = "replicate", value = "luminescence",
                              normalized_r1_ch1:normalized_r3_ch1)
samples_norm <- samples_norm %>% mutate(BRCA1 = case_when(plate == 1 ~ "Mutant",
                                                plate == 2 ~ "Restored"))
hn <- ggplot(samples_norm, aes(x=luminescence, fill=BRCA1, color=BRCA1)) +
  geom_histogram(position="identity", bins = 40, alpha=0.5) +
  ylab("Count") + xlab("Raw Luminescence")
hn
