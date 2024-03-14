### Fisher's test for diff spliced genes in mnat

if(!require(ggstatsplot)){
  install.packages("ggstatsplot")
}

library(ggstatsplot)
library(cowplot)

Genes_with_RA = 52

RA_as_01 = 13

RA_as_01/Genes_with_RA # [1] 0.25

Genes_without_RA = 16085

nonRA_as_01 = 1361

nonRA_as_01/Genes_without_RA # [1] 0.08461299

a = matrix(c(Genes_with_RA, RA_as_01, Genes_without_RA, nonRA_as_01), nrow = 2)

res = fisher.test(a)

df = data.frame(c(rep('Rearranged', Genes_with_RA), rep('Non-rearranged', Genes_without_RA)), 
           c(rep('AS', RA_as_01), rep('non-AS', Genes_with_RA - RA_as_01), 
             rep('AS', nonRA_as_01), rep('non-AS', Genes_without_RA - nonRA_as_01)))

names(df) = c('RA', 'AS')

table(df$RA)
table(df$AS)

res_plot = ggbarstats(
  df, AS, RA,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(res$p.value < 0.001, "< 0.001", round(res$p.value, 3))
  )
)

save_plot("FisherTest_01_mnat.pdf",
          res_plot)

### RA in the origin of true workers

Genes_with_RA_advanced = 149
RA_as_01 = 33

Genes_without_RA_advanced = 15988
nonRA_as_01 = 1341

ct_advanced = matrix(c(Genes_with_RA_advanced, RA_as_01, Genes_without_RA_advanced, nonRA_as_01), nrow = 2)

res_advanced = fisher.test(ct_advanced)

res_advanced

df_advanced = data.frame(c(rep('Rearranged', Genes_with_RA_advanced), rep('Non-rearranged', Genes_without_RA_advanced)), 
                c(rep('AS', RA_as_01), rep('non-AS', Genes_with_RA_advanced - RA_as_01), 
                  rep('AS', nonRA_as_01), rep('non-AS', Genes_without_RA_advanced - nonRA_as_01)))

names(df_advanced) = c('RA', 'AS')

table(df_advanced$RA)
table(df_advanced$AS)


res_advanced_plot = ggbarstats(
  df_advanced, AS, RA,
  results.subtitle = FALSE,
  subtitle = paste0(
    "Fisher's exact test", ", p-value = ",
    ifelse(res_advanced$p.value < 0.001, "< 0.001", round(res_advanced$p.value, 3))
  )
)

save_plot("FisherTest_01_advanced_mnat.pdf",
          res_advanced_plot)






calculate_fisher_AS = function(Genes_with_RA, Genes_without_RA, RA_as, nonRA_as, output){
  ct = matrix(c(Genes_with_RA, RA_as, Genes_without_RA, nonRA_as), nrow = 2)
  
  res = fisher.test(ct)
  
  res
  
  df = data.frame(c(rep('Rearranged', Genes_with_RA), rep('Non-rearranged', Genes_without_RA)), 
                           c(rep('AS', RA_as), rep('non-AS', Genes_with_RA - RA_as), 
                             rep('AS', nonRA_as), rep('non-AS', Genes_without_RA - nonRA_as)))
  
  names(df) = c('RA', 'AS')
  
  table(df$RA)
  table(df$AS)
  
  
  res_plot = ggbarstats(
    df, AS, RA,
    results.subtitle = FALSE,
    subtitle = paste0(
      "Fisher's exact test", ", p-value = ",
      ifelse(res$p.value < 0.001, "< 0.001", round(res$p.value, 3))
    )
  )
  
  save_plot(output, res_plot)
}

# kings as_05 termite origin
calculate_fisher_AS(52, 16085, 13, 1538, "king_05_fisher.pdf")

# T0 as_05 termite origin
calculate_fisher_AS(52, 16085, 9, 1330, "T0_05_fisher.pdf")

# T4 as_05 termite origin
calculate_fisher_AS(52, 16085, 3, 233, "T4_05_fisher.pdf")

# T4 as_05 termite origin
calculate_fisher_AS(52, 16085, 6, 364, "worker_05_fisher.pdf")

# kings as_05 both origins
calculate_fisher_AS(201, 15936, 41, 1510, "king_05_fisher_bothRA.pdf")

# T0 as_05 both origins
calculate_fisher_AS(201, 15936, 35, 1304, "T0_05_fisher_bothRA.pdf")

# T4 as_05 both origins
calculate_fisher_AS(201, 15936, 6, 230, "T4_05_fisher_bothRA.pdf")

# worker as_05 both origins
calculate_fisher_AS(201, 15936, 13, 357, "worker_05_fisher_bothRA.pdf")
