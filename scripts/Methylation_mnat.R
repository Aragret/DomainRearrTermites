if(!require(dplyr)){
  install.packages("dplyr")
}

if(!require(tidyr)){
  install.packages("tidyr")
}

if(!require(ggplot2)){
  install.packages("ggplot2")
}

if(!require(ggpubr)){
  install.packages("ggpubr")
}


library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)

#### INPUT ####

Methylation_table = 'data/mnat_methylation/mnat_methylation'
Rearranged_table = 'analysis/domrates/termites_arrangements_exact.out'
Advanced_RA_table = 'analysis/domrates/advanced_sociality_arrangements_exact.out'
log_file = 'mnat_methylation_comparison.out'
Annotation_file = 'Mnat_pep_noPG_loI.dom' # from the Rearranged_table
Output_Figure = 'mnat_methylation_comparison.pdf'

######################

# read tables
Methylation = read.table(Methylation_table, sep='\t', header = TRUE)
rearr_genes = read.table(Rearranged_table, sep='\t')
advanced_rearr_genes = read.table(Advanced_RA_table, sep='\t')

all_RA = rbind(rearr_genes, advanced_rearr_genes)
# rename columns
names(all_RA) = c('AnnotationFile', 'Node', 'Exact', 'Type', 'Arrangement', 'Genes')

# modify RA table to one line - one gene
sp_rearr_genes = all_RA %>% 
  filter(AnnotationFile == Annotation_file) %>%
  mutate(Gene = strsplit(Genes, ","),
         Rearr = 1) %>%
  unnest(Gene) %>%
  select(Gene, Rearr)

summary(Methylation$Mean)

Mean_methylation = Methylation %>%
  group_by(Gene) %>%
  # filter(Feature == 'cds') %>%
  mutate(Same_NumberSites = NumberSites - mean(NumberSites)) %>%
  # filter(Same_NumberSites == 0) %>%
  summarise(Gene_Mean = mean(Mean)) %>%
  mutate(Gene_Mean_log = log(Gene_Mean))

df = left_join(Mean_methylation, sp_rearr_genes) %>%
  mutate(RA = as.factor(replace_na(Rearr, 0)))

table(Methylation$Feature)

wilcox.test(df[df$RA == 0,]$Gene_Mean, df[df$RA == 1,]$Gene_Mean)

summary(df[df$RA == 1,]$Gene_Mean)
summary(df[df$RA == 0,]$Gene_Mean)

sink(log_file)

print(Sys.Date())
print('Log of the Methylation_mnat.R script, comparing methylation for genes rearranged in the origin of termites and true workers')

print(c('N of rearranged:', nrow(df[df$RA == 1,]))) # N for RA

print(c('N of non-rearranged:', nrow(df[df$RA == 0,]))) # N for non-RA

print(c('Stats for non-rearranged methylation (mean and var):', mean(df[df$RA == 0,]$Gene_Mean), 
        var(df[df$RA == 0,]$Gene_Mean)))
print(c('Stats for rearranged methylation (mean and var):', mean(df[df$RA == 1,]$Gene_Mean), 
        var(df[df$RA == 1,]$Gene_Mean)))

print('wilcox.test one-sided')
print(wilcox.test(df[df$RA == 0,]$Gene_Mean, df[df$RA == 1,]$Gene_Mean, alternative = 'l'))

sink()

pdf(Output_Figure)

gghistogram(df, x = "Gene_Mean") + 
  facet_grid(rows = vars(RA), scales = 'free_y', labeller = label_both) +
  labs(x = 'Mean methylation')

gghistogram(df, x = "Gene_Mean_log") + 
  facet_grid(rows = vars(RA), scales = 'free_y', labeller = label_both) +
  labs(x = 'Mean methylation, log')

dev.off()

ggdensity(df, x = "Gene_Mean", y = "density", color = "RA", add = "mean") + 
  labs(x = 'Mean methylation')

# write.csv(df, "MnatMethylation.csv", row.names = FALSE, quote = FALSE)
