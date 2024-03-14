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

#### CHANGE INPUT ####

CpGoe_table = 'dpun_CpGoe_alldom.tsv'
Rearranged_table = 'analysis/domrates/termites_ancestral_arrangements'
RA_advanced = 'analysis/domrates/advanced_sociality_ancestral_arrangements'
log_file = 'dpun_CpGoe_comparison.out'
Annotation_file = 'Dpun_OGSprelim3_pep.dom' # from the Rearranged_table
Output_Figure = 'dpun_CpGoe_comparison.pdf'

######################

CpG_comparison = function(CpGoe_table, Rearranged_table, log_file, Annotation_file, Output_Figure){
  
  # read tables
  CpGoe = read.table(CpGoe_table, sep='\t', header = TRUE)
  rearr_genes = read.table(Rearranged_table, sep='\t')
  
  # rename columns
  names(rearr_genes) = c('AnnotationFile', 'Node', 'Exact', 'Type', 'Arrangement', 'Genes')
  
  # modify RA table to one line - one domain RA
  sp_rearr_genes = rearr_genes %>% 
    filter(AnnotationFile == Annotation_file & 
             Exact == 'exact solution' &
             Type != 'fusion') %>%
    mutate(GeneID = strsplit(Genes, ",")) %>%
    unnest(GeneID) %>% 
    mutate(Domain = strsplit(Arrangement, "-")) %>%
    unnest(Domain)
  
  # extract geneID and domain names in CpGoe table
  
  CpGoe = CpGoe %>%
    separate(Name, c('Name', 'Domain'), ' ') %>%
    separate(Name, c('GeneID1', 'GeneID2', 'GeneID3', NA), '_') %>% ### PLEASE CHANGE THE SEPARATE RULES depending on the input file
    unite('GeneID', GeneID1:GeneID3)
  
  # merge the tables
  df = full_join(CpGoe, sp_rearr_genes[, c('GeneID', 'Domain', 'Arrangement', 'Type')]) %>%
    mutate(Rearranged = case_when(is.na(Arrangement) ~ 0,
                                  !is.na(Arrangement) ~ 1))
  
  # calculate the weighted mean for each arrangement
  mean_CpGoe_RA = df %>%
    filter(Rearranged == 1) %>%
    group_by(GeneID, Arrangement) %>%
    summarise(CpG_mean = weighted.mean(CpGoe, CG)) %>%
    mutate(Rearranged = 1)
  
  mean_CpGoe_nonRA = df %>%
    filter(Rearranged == 0) %>%
    group_by(GeneID) %>%
    summarise(CpG_mean = weighted.mean(CpGoe, CG)) %>%
    mutate(Rearranged = 0,
           Arrangement = NA)
  
  mean_CpGoe = bind_rows(mean_CpGoe_nonRA, mean_CpGoe_RA) %>%
    mutate(Rearranged = as.factor(Rearranged),
           CpG_mean_log = log(CpG_mean))
  
  # add important info to the log file
  sink(log_file)
  
  print(Sys.Date())
  print(c('Log of the CpGoe_comparison.R script for', Annotation_file))
  
  print('Summary of C+G')
  print(summary(CpGoe$CG)) # number of frC + frG for the weighted average
  
  print('Summary of domain lengths')
  print(summary(CpGoe$Length)) # number of frC + frG for the weighted average
  
  print(c('N of rearranged:', nrow(mean_CpGoe_RA))) # N for RA
  
  print(c('N of non-rearranged:', nrow(mean_CpGoe_nonRA))) # N for non-RA
  
  print(c('Stats for non-rearranged CpGoe (mean and var):', mean(mean_CpGoe_nonRA$CpG_mean), 
          var(mean_CpGoe_nonRA$CpG_mean)))
  print(c('Stats for rearranged CpGoe (mean and var):', mean(mean_CpGoe_RA$CpG_mean), 
          var(mean_CpGoe_RA$CpG_mean)))
  
  print('wilcox.test(l) one-sided')
  print(wilcox.test(mean_CpGoe_RA$CpG_mean, mean_CpGoe_nonRA$CpG_mean, alternative = 'l'))
  
  sink()
  
  pdf(Output_Figure)
  
  print(gghistogram(mean_CpGoe, x = "CpG_mean_log") + 
          facet_grid(rows = vars(Rearranged), scales = 'free_y', labeller = label_both))
  
  dev.off()
  
}

CpG_comparison_allRA = function(CpGoe_table, Rearranged_table, RA_advanced, log_file, Annotation_file, Output_Figure){
  
  # read tables
  CpGoe = read.table(CpGoe_table, sep='\t', header = TRUE)
  rearr_genes = read.table(Rearranged_table, sep='\t')
  RA_advanced = read.table(RA_advanced, sep='\t')
  
  # rename columns
  names(rearr_genes) = c('AnnotationFile', 'Node', 'Exact', 'Type', 'Arrangement', 'Genes')
  names(RA_advanced) = c('AnnotationFile', 'Node', 'Exact', 'Type', 'Arrangement', 'Genes')
  
  # modify RA table to one line - one domain RA
  sp_rearr_genes = rearr_genes %>% 
    filter(AnnotationFile == Annotation_file & 
             Exact == 'exact solution' &
             Type != 'fusion') %>%
    mutate(GeneID = strsplit(Genes, ",")) %>%
    unnest(GeneID) %>% 
    mutate(Domain = strsplit(Arrangement, "-")) %>%
    unnest(Domain)
  
  sp_RA_advanced = RA_advanced %>% 
    filter(AnnotationFile == Annotation_file & 
             Exact == 'exact solution' &
             Type != 'fusion') %>%
    mutate(GeneID = strsplit(Genes, ",")) %>%
    unnest(GeneID) %>% 
    mutate(Domain = strsplit(Arrangement, "-")) %>%
    unnest(Domain)
  
  all_RA = bind_rows(sp_rearr_genes, sp_RA_advanced)
  
  # extract geneID and domain names in CpGoe table
  
  CpGoe = CpGoe %>%
    separate(Name, c('Name', 'Domain'), ' ') %>%
    separate(Name, c('GeneID', NA), '_') ### PLEASE CHANGE THE SEPARATE RULES depending on the input file
  # unite('GeneID', GeneID1:GeneID3)
  
  # merge the tables
  df = full_join(CpGoe, all_RA[, c('GeneID', 'Domain', 'Arrangement', 'Type')]) %>%
    mutate(Rearranged = case_when(is.na(Arrangement) ~ 0,
                                  !is.na(Arrangement) ~ 1))
  
  # calculate the weighted mean for each arrangement
  mean_CpGoe_RA = df %>%
    filter(Rearranged == 1) %>%
    group_by(GeneID, Arrangement) %>%
    summarise(CpG_mean = weighted.mean(CpGoe, CG)) %>%
    mutate(Rearranged = 1)
  
  mean_CpGoe_nonRA = df %>%
    filter(Rearranged == 0) %>%
    group_by(GeneID) %>%
    summarise(CpG_mean = weighted.mean(CpGoe, CG)) %>%
    mutate(Rearranged = 0,
           Arrangement = NA)
  
  mean_CpGoe = bind_rows(mean_CpGoe_nonRA, mean_CpGoe_RA) %>%
    mutate(Rearranged = as.factor(Rearranged),
           CpG_mean_log = log(CpG_mean))
  
  # add important info to the log file
  sink(log_file)
  
  print(Sys.Date())
  print(c('Log of the CpGoe_comparison.R script for', Annotation_file))
  
  print('Summary of C+G')
  print(summary(CpGoe$CG)) # number of frC + frG for the weighted average
  
  print('Summary of domain lengths')
  print(summary(CpGoe$Length)) # number of frC + frG for the weighted average
  
  print(c('N of rearranged:', nrow(mean_CpGoe_RA))) # N for RA
  
  print(c('N of non-rearranged:', nrow(mean_CpGoe_nonRA))) # N for non-RA
  
  print(c('Stats for non-rearranged CpGoe (mean and var):', mean(mean_CpGoe_nonRA$CpG_mean), 
          var(mean_CpGoe_nonRA$CpG_mean)))
  print(c('Stats for rearranged CpGoe (mean and var):', mean(mean_CpGoe_RA$CpG_mean), 
          var(mean_CpGoe_RA$CpG_mean)))
  
  print('wilcox.test(l) one-sided')
  print(wilcox.test(mean_CpGoe_RA$CpG_mean, mean_CpGoe_nonRA$CpG_mean, alternative = 'l'))
  
  sink()
  
  pdf(Output_Figure)
  
  print(gghistogram(mean_CpGoe, x = "CpG_mean_log") + 
          facet_grid(rows = vars(Rearranged), scales = 'free_y', labeller = label_both))
  
  dev.off()
  
}

CpG_comparison(CpGoe_table, Rearranged_table, log_file, Annotation_file, Output_Figure)

CpG_comparison_allRA(CpGoe_table, Rearranged_table, RA_advanced,
                     log_file = 'pame_allRA_CpGoe_comparison.out',
                     Annotation_file = Annotation_file,
                     Output_Figure = 'pame_allRA_CpGoe_comparison.pdf')

