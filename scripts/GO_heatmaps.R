# GO enrichment of genes with domain rearrangements

library("topGO")
library("RColorBrewer")
library("tagcloud")
library(ViSEAGO)
library(GO.db)
library(GOSemSim)
library(ggdendro)
library(grid)
library(ggplot2)

# load the GO result tables for all species for visualisation
allRes_BP_cfor = read.table("GO_analysis/Cfor_advanced_BP_table",
                            header = TRUE, sep = "\t")

allRes_BP_csec$Species = "csec"
allRes_BP_znev$Species = "znev"
allRes_BP_rspe_advanced$Species = "rspe"
allRes_BP_mnat_advanced$Species = "mnat"
allRes_bger$Species = "bger"
allRes_pame$Species = "pame"
allRes_BP_cfor$Species = "cfor"
allRes_dpun$Species = "dpun"



df = rbind(allRes_BP_csec, allRes_BP_znev, allRes_BP_rspe_advanced, allRes_BP_mnat_advanced, allRes_bger, 
           allRes_pame, allRes_BP_cfor[, -7], allRes_dpun)

df$Species = factor(df$Species, levels = c("bger", "dpun", "pame", "znev", "csec", "cfor", "rspe", "mnat"))

df[df$pvalue == "< 1e-30",]$pvalue = 1e-30

ggplot(df, aes(Species, Term, fill = as.numeric(pvalue))) +
  geom_tile() + 
  scale_fill_gradient(low="blue", high="white")


### adding semantic similarity

dmGO <- godata('org.Dm.eg.db', ont="BP")

# a = mgoSim(allRes_BP_csec$GO.ID, allRes_BP_znev$GO.ID, semData=hsGO, measure="Wang", combine=NULL)
# 
# a = mgoSim(allRes_BP_csec$GO.ID, allRes_BP_znev$GO.ID, semData=hsGO, measure="Wang", combine=NULL)
# 
# b = termSim(df$Term, df$Term, semData=hsGO)

# clustering GO terms

SemSim_m = mgoSim(unique(df$GO.ID), unique(df$GO.ID), semData = dmGO, combine = NULL)


SemSim_dendro <- as.dendrogram(hclust(d = dist(x = SemSim_m)))

dendro_plot <- ggdendrogram(data = SemSim_dendro, rotate = TRUE, labels = FALSE,
                            leaf_labels = FALSE)

dendro_plot + theme(axis.text.y = element_text(size = 6))

df$GO.ID <- factor(x = df$GO.ID,
                               levels = colnames(SemSim_m)[order.dendrogram(SemSim_dendro)], 
                               ordered = TRUE)

df[df$Term == "protein modification by small protein co...",]

df[81,]$Term = "protein modification by small protein conjugation"
df[83,]$Term = "protein modification by small protein conjugation or removal"

colnames(SemSim_m)[order.dendrogram(SemSim_dendro)]
length(colnames(SemSim_m)[order.dendrogram(SemSim_dendro)])

# ordered_terms = c("DNA repair", "RNA processing", "DNA metabolic process", "nucleic acid metabolic process",
#                   "RNA metabolic process", "protein phosphorylation", "protein ubiquitination", "protein modification by small protein conjugation",
#                   "protein modification by small protein conjugation or removal", "gene expression",
#                   "macromolecule modification", "macromolecule metabolic process", "organonitrogen compound metabolic proces...",
#                   "protein metabolic process", "protein modification process", "monosaccharide metabolic process",
#                   "hexose metabolic process", "mannose metabolic process", "organic substance metabolic process",
#                   "metabolic process", "primary metabolic process", "nitrogen compound metabolic process",
#                   "dephosphorylation", "phosphorylation", "phosphate-containing compound metabolic ...",
#                   "phosphorus metabolic process", "cellular metabolic process", "regulation of cellular process",
#                   "cell communication", "cellular process", "establishment of localization in cell",
#                   "intracellular transport", "cellular macromolecule localization", "cellular localization",
#                   "protein-containing complex assembly", "cell-cell adhesion via plasma-membrane a...",
#                   "homophilic cell adhesion via plasma memb...", "cell adhesion", "cell-cell adhesion",
#                   "multicellular organismal process", "signaling", "positive regulation of metabolic process",
#                   "regulation of metabolic process", "regulation of multicellular organismal p...",
#                   "biological regulation", "regulation of biological process", "positive regulation of biological proces...",
#                   "regulation of ARF protein signal transdu...", "ARF protein signal transduction",
#                   "G protein-coupled receptor signaling pat...", "intracellular signal transduction",
#                   "signal transduction", "response to stress", "response to stimulus", "cellular response to stimulus",
#                   "cellular response to stress", "cellular response to DNA damage stimulus",
#                   "establishment of protein localization", "protein localization", "macromolecule localization",
#                   "intracellular protein transport", "protein transport", "nitrogen compound transport",
#                   "carbohydrate transport", "organic substance transport", "cell development",
#                   "cell differentiation", "cellular developmental process", "neuron development",
#                   "neuron differentiation", "neurogenesis", "generation of neurons", 
#                   "multicellular organism development", "nervous system development", "system development")

ordered_terms = c("DNA repair", "transposition, DNA-mediated", "DNA recombination", "DNA integration",
                  "RNA processing", "RNA metabolic process", "nucleic acid metabolic process", 
                  "protein glycosylation", "glycoprotein biosynthetic process", "glycoprotein metabolic process",
                  "protein phosphorylation", "protein modification by small protein conjugation or removal", 
                  "protein ubiquitination", "protein modification by small protein conjugation", "macromolecule glycosylation",
                  "gene expression", "macromolecule modification", "macromolecule metabolic process", 
                  "organonitrogen compound metabolic proces...", "protein metabolic process", "protein modification process", 
                  "monosaccharide metabolic process", "mannose metabolic process", "hexose metabolic process",
                  "glycosylation", "primary metabolic process", "nitrogen compound metabolic process",
                  "organic substance metabolic process", "metabolic process", "intracellular protein transport", 
                  "protein transport", "organic substance transport", "carbohydrate transport", 
                  "nitrogen compound transport", "ion transport", "transport", "cellular macromolecule localization", 
                  "cellular localization", "intracellular transport", "establishment of localization in cell",
                  "transmembrane transport", "Golgi vesicle transport", "establishment of localization",
                  "localization", "establishment of protein localization", "protein localization", 
                  "macromolecule localization", "cellular process", "transposition", "cell communication",
                  "dephosphorylation", "phosphorylation", "phosphate-containing compound metabolic ...",
                  "phosphorus metabolic process", "cellular metabolic process", "positive regulation of metabolic process",
                  "regulation of metabolic process", "regulation of cellular process",
                  "biological regulation", "regulation of biological process", "positive regulation of biological proces...",
                  "cell projection organization", "plasma membrane bounded cell projection ...",
                  "homophilic cell adhesion via plasma memb...", "cell-cell adhesion via plasma-membrane a...",
                  "cell adhesion", "cell-cell adhesion", "regulation of ARF protein signal transdu...", "ARF protein signal transduction",
                  "G protein-coupled receptor signaling pat...", "signal transduction", "intracellular signal transduction",
                  "cellular response to stimulus", "cellular response to DNA damage stimulus", "cellular response to stress", 
                  "signaling", "response to stress", "response to stimulus")

                  
                  # "cell communication", 
                  # "protein-containing complex assembly", 
                  # "multicellular organismal process", 
                  # "regulation of multicellular organismal p...",
                  # 
                  # 
                  # 
                  # 
                  # 
                  # "cell development",
                  # "cell differentiation", "cellular developmental process", "neuron development",
                  # "neuron differentiation", "neurogenesis", "generation of neurons", 
                  # "multicellular organism development", "nervous system development", "system development")

df$Term <- factor(x = df$Term,
                   levels = rev(ordered_terms), 
                   ordered = TRUE)

heatmap_plot = ggplot(df, aes(Species, Term, fill = as.numeric(pvalue))) +
  geom_tile() + 
  scale_fill_gradient(low="blue", high="white") +
  theme(axis.text.y = element_text(size = 8),
        legend.position = "top")


pdf("TermiteWorkersOriginsGenesRoaches_heatmap.pdf",
    width = 8)

grid.newpage()
print(heatmap_plot, 
      vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(dendro_plot, 
      vp = viewport(x = 0.90, y = 0.43, width = 0.2, height = 0.97))

dev.off()
