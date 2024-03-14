#!/usr/bin/env Rscript

# Modified by Alina to add the visualization

suppressMessages(library("topGO"))
suppressMessages(library("tagcloud"))
suppressMessages(library("RColorBrewer"))

args = commandArgs(trailingOnly=TRUE)

# load the GO universe
GOmap <- readMappings(file = args[1])

# args
# 1: universe
# 2: genes of interest
# 3: algorithm
# 4: output

# make list of gene names in the universe
genes_in_universe <- names(GOmap)

# loading and preparing of the genes of interest
subset <- scan(file=args[2], what=character())
genes_of_interest <- factor(as.integer(genes_in_universe %in% subset))
names(genes_of_interest) <- genes_in_universe

# performing the enrichment analysis for biological processes
GOdata_BP <- new("topGOdata", description="GO-Analysis", ontology="BP", allGenes=genes_of_interest,  annot = annFUN.gene2GO, gene2GO = GOmap)
GOdata_MF <- new("topGOdata", description="GO-Analysis", ontology="MF", allGenes=genes_of_interest,  annot = annFUN.gene2GO, gene2GO = GOmap)

# run a significance test
resultParentchild_BP <- runTest(GOdata_BP, algorithm=args[3], statistic="fisher") 
resultParentchild_MF <- runTest(GOdata_MF, algorithm=args[3], statistic="fisher") 

mf_out = paste(args[4], "_MF.txt", sep = "")
bp_out = paste(args[4], "_BP.txt", sep = "")

# filter results for GO terms with a p-value <= 0.05
summary_BP <- summary(attributes(resultParentchild_BP)$score <= 0.05)
numsignif_BP <- as.integer(summary_BP[[3]])

if(numsignif_BP > 0)
{
	results_BP <- GenTable(GOdata_BP, fisher = resultParentchild_BP, topNodes = numsignif_BP)
} else {
	results_BP <- NULL
}

write.csv(results_BP, file=bp_out)

if(length(results_BP$Term) > 1){
  for(i in 1:length(results_BP$fisher)){
    if(results_BP$fisher[i] == "< 1e-30"){
      results_BP$fisher[i] <- 1e-30
    }
  }

  if(length(results_BP$Term) > 20){
    results_BP <- results_BP[1:20,]
  }


  tag_BP <- paste0(args[4], "_BP_tagcloud.pdf")

  pdf(file = tag_BP,width=30, height = 5)
  colors <- colorRampPalette( brewer.pal( 12, "Paired" ) )( length(results_BP$Term) )
  tagcloud(strmultline(results_BP$Term), weights = -log(as.numeric(results_BP$fisher)),col=colors)
  dev.off()
}


summary_MF <- summary(attributes(resultParentchild_MF)$score <= 0.05)
numsignif_MF <- as.integer(summary_MF[[3]])

if(numsignif_MF > 0)
{
	results_MF <- GenTable(GOdata_MF, fisher = resultParentchild_MF, topNodes = numsignif_MF)
} else {
	results_MF <- NULL
}
	

# write results to file in simple CSV format
write.csv(results_MF, file=mf_out)

if(length(results_MF$Term) > 1){
  for(i in 1:length(results_MF$fisher)){
    if(results_MF$fisher[i] == "< 1e-30"){
      results_MF$fisher[i] <- 1e-30
    }
  }

  if(length(results_MF$Term) > 20){
    results_MF <- results_MF[1:20,]
  }

  tag_MF <- paste0(args[4], "_MF_tagcloud.pdf")

  pdf(file = tag_MF,width=30, height = 5)
  colors <- colorRampPalette( brewer.pal( 12, "Paired" ) )( length(results_MF$Term) )
  tagcloud(strmultline(results_MF$Term), weights = -log(as.numeric(results_MF$fisher)),col=colors)
  dev.off()

}

