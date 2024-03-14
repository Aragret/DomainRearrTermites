library(dplyr)
library(ggstatsplot)

###### INPUT
# change bger to other species to complete the analysis on the rest of the species
TEcounts_table = 'data/TE_proximity/bger_summary_mod.csv'
RA_genes_table = 'analysis/domrates/cockroaches_RA_genes_mod'
output_pdf = 'bger_out.pdf'


calculate_TE_RA = function(TEcounts_table, RA_genes_table, output_pdf){

  # read the input tables
  TEcounts = read.csv(TEcounts_table)
  RA_genes = read.table(RA_genes_table)
  
  names(RA_genes) = 'gene_id'
  
  RA_genes$RA = 1
  
  TEcounts_all = TEcounts %>%
    select(gene_id, gene_all.TE., promotor_all.TE., X10kb.flanks_all.TE.) %>%
    left_join(RA_genes) %>%
    replace(is.na(.), 0)
  
  # wilcox.test(TEcounts_all[TEcounts_all$RA == 0,]$gene_all.TE., TEcounts_all[TEcounts_all$RA == 1,]$gene_all.TE.,
  #             alternative = "greater")
  # 
  # wilcox.test(TEcounts_all[TEcounts_all$RA == 0,]$promotor_all.TE., 
  #             TEcounts_all[TEcounts_all$RA == 1,]$promotor_all.TE.,
  #             alternative = "greater")
  # 
  # wilcox.test(TEcounts_all[TEcounts_all$RA == 0,]$X10kb.flanks_all.TE., 
  #             TEcounts_all[TEcounts_all$RA == 1,]$X10kb.flanks_all.TE.,
  #             alternative = "greater")
  
  set.seed(123)
  
  pdf(output_pdf)
  
  print(ggbetweenstats(
    data = TEcounts_all,
    x = RA,
    y = gene_all.TE.,
    type = "nonparametric",
    title = "Distribution of TE counts across non-rearrangend and rearranged genes"
  ))
  
  print(ggbetweenstats(
    data = TEcounts_all,
    x = RA,
    y = promotor_all.TE.,
    type = "nonparametric",
    title = "Distribution of TE counts across non-rearrangend and rearranged genes"
  ))
  
  print(ggbetweenstats(
    data = TEcounts_all,
    x = RA,
    y = X10kb.flanks_all.TE.,
    type = "nonparametric",
    title = "Distribution of TE counts across non-rearrangend and rearranged genes"
  ))
  
  dev.off()
}

calculate_TE_RA(TEcounts_table, RA_genes_table, output_pdf)

