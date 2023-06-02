setwd("C:/Users/misak/Desktop/ICR/siRNA")
library(dplyr)

all_six_z <- read.csv("all_six_zscores.txt", sep="\t", header=TRUE)
z_diff <- read.csv("zscores_difference.txt", sep="\t", header=TRUE)

# Welch's t-test for all genes (pairwise comparison)
library(matrixTests)
cols <- as.vector(colnames(all_six_z))
mutcols <- cols[4:9] # storing column names for mutant z scores
rescols <- cols[10:15] # storing column names for restored z scores
t_result <- row_t_welch(all_six_z[,mutcols], all_six_z[,rescols]) %>% 
  mutate(gene = all_six_z$gene) # dataframe with t test result for each row containing just the z scores+ gene names

# Export dataframe with six z scores and p values
pval <- t_result %>% select(pvalue, gene)
z_diff <- left_join(z_diff, pval, by = "gene")
write.table(
  z_diff,
  "z_score_difference_pval.txt",
  col.names=TRUE,
  sep="\t",
  quote=FALSE,
  row.names=FALSE
)


# Group genes based on p-value range
pval <- t_result %>% select(pvalue, gene)
p0.01 <- pval %>% filter(pvalue < 0.01) %>% 
  mutate(pstatus = "p<0.01")
p0.05 <- pval %>% filter(pvalue < 0.05 & pvalue >= 0.01) %>% 
  mutate(pstatus = "p<0.05") 
p0.1 <- pval %>% filter(pvalue < 0.1 & pvalue >= 0.05) %>% 
  mutate(pstatus = "p<0.1")
pNS <- pval %>% filter(pvalue >= 0.1) %>% 
  mutate(pstatus = "N.S.")
all_pstatus <- bind_rows(p0.01, p0.05, p0.1, pNS)

z_diff <- left_join(z_diff, all_pstatus, by = "gene")


# Change pstatus factor levels so that p<0.01 will have darkest color (e.g. bottom of the legend)
levels(z_diff$pstatus)
z_diff$pstatus <- factor(z_diff$pstatus,
                         levels = c("N.S.", "p<0.1", "p<0.05", "p<0.01"))

# Waterfall with bars ordered in increasing mut_res_diff


p <- ggplot(z_diff, aes(x = reorder(gene, +mut_res_diff), 
                        y = mut_res_diff,
                        fill = pstatus)) +
  theme_bw() +
  geom_bar(stat="identity", width = 1) +
  scale_fill_brewer(palette = "Greens") +
  theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))
p

# new_z_diff has gene names only when p<0.01

new_z_diff <- z_diff %>% mutate(gene_name = case_when(
  pstatus == "N.S." ~ " ",
  pstatus == "p<0.1" | pstatus == "p<0.05" | pstatus == "p<0.01" ~ gene
))
new_z_diff <- new_z_diff %>% arrange(mut_res_diff)
goi <- as.vector(new_z_diff$gene_name)

p <- ggplot(new_z_diff, aes(x = reorder(gene, +mut_res_diff), 
                        y = mut_res_diff,
                        fill = pstatus)) +
  theme_bw(base_size = 15) +
  theme(legend.position = c(0.5, 0.1),
        legend.direction="horizontal",
        legend.title=element_blank(), 
        axis.text.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size=.1),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank()) +
  geom_bar(stat="identity", width = 1) +
  scale_fill_brewer(palette = "Greens") +
  ylab("z-score difference") + xlab("") +
  geom_col() +
  geom_text(size = 3.5, aes(y = 0, label = gene_name, 
                hjust = ifelse(mut_res_diff < 0, 0, 1), # ifelse(test, yes, no)
                vjust = 0.5,
                angle = 90))
p
# hjust, vjust reference: https://stackoverflow.com/questions/7263849/what-do-hjust-and-vjust-do-when-making-a-plot-using-ggplot


             