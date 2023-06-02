library(dplyr)
library(ggplot2)
library(ggrepel)
setwd("C:/Users/misak/Desktop/ICR/siRNA")

# load dataframe with six z-scores and p-value
z_orip <- read.csv("z_score_difference_pval.txt", sep="\t", header=TRUE)

# calculate adjusted p-value
ori_pval <- as.vector(z_orip$pvalue)
adj_pval <- p.adjust(ori_pval, method = 'BH') #Benjamini-Hochberg method
z_adjp <- z_orip
z_adjp$adj_pval <- adj_pval

# volcano plot


p1 <- ggplot(z_orip, aes(mut_res_diff, -log(adj_pval,10))) +   
  geom_point(size = 1) +
  xlab('Z-score difference') + 
  ylab(expression("-log"[10]*"p-value")) +
  xlim(-4,4)
p1

# add color to statistically significant (p < 0.05)

data_color <- z_adjp %>% 
  mutate(
    Sensitivity = case_when(mut_res_diff >= 1 & adj_pval <= 0.05 ~ "Restored sensitive",
                           mut_res_diff <= -1 & adj_pval <= 0.05 ~ "Mutant sensitive",
                           TRUE ~ "N.S.")
  )
p2 <- ggplot(data_color, aes(mut_res_diff, -log(adj_pval,10))) +
  theme_bw() +
  geom_point(aes(color = Sensitivity), size = 1.5) +
  xlab('Z-score difference') + 
  ylab(expression("-log"[10]*"p-value")) +
  xlim(-4, 4) +
  scale_color_manual(values = c("dodgerblue3", "gray50", "firebrick3")) + # legends are ordered alphabetically
  guides(colour = guide_legend(override.aes = list(size=2))) # without override, legend dots will be the same size as plotted dots
p2

# add gene label to statistically significant hits

top <- 10
top_genes <- bind_rows( # make new dataset with just top 10 differentially over/under-expressed genes
  data_color %>% 
    filter(Sensitivity == 'Restored sensitive') %>% 
    arrange(adj_pval, desc(abs(mut_res_diff))) %>% 
    head(top),
  data_color %>% 
    filter(Sensitivity == 'Mutant sensitive') %>% 
    arrange(adj_pval, desc(abs(mut_res_diff))) %>% 
    head(top)
)


p3 <-  p2 +
  geom_label_repel(data = top_genes,
                   mapping = aes(mut_res_diff, -log(adj_pval,10), label = gene),
                   size = 2) +
  geom_point(data=top_genes,
             aes(mut_res_diff, -log(adj_pval,10)), 
             size = 3, shape = 21) +
  geom_hline(yintercept = -log(0.05,10), linetype="dashed", 
             color = "grey10", linewidth=0.5) +
  geom_vline(xintercept = -1, linetype="dashed", 
             color = "grey10", linewidth=0.5) +
  geom_vline(xintercept = 1, linetype="dashed", 
             color = "grey10", linewidth=0.5)
p3

