## Overview

This repository contains data and code used in the article "Major changes in domain arrangements are associated with the evolution of termites" (BioRxiv version: https://doi.org/10.1101/2023.05.15.540413, accepted in JEB).

## Repository structure

* data/proteomes contains publically available proteomes of Blattodea (for references please check out the Methods section of the paper above)

* data/pfam contains the annotated domains for the proteomes by using pfam_scan.pl version 1.6 database version v33.1

* scripts/ contains code used in the analyses

* analyis/ contains intermediate files for each analysis

## Workflow

### Filtering

The downloaded proteomes were filtered to keep the longest isoforms and get rid of pseudogenes using scripts from DW-Helper: https://domain-world.zivgitlabpages.uni-muenster.de/dw-helper/

``` isoformCleaner -i proteome.fa -o proteome_loI.fa ```

``` seqCheck -i proteome_loI.fa -o proteome_loI_noPG.fa --fix-and-remove ```

The completeness of the proteomes was checked with the DOGMA webserver: https://domainworld-services.uni-muenster.de/dogma/ (insects set, Pfam 35 version). To check whether different groups of interest show significant differences in the proteome completeness, the script ``` scripts/anova_DOGMA.R ``` was used.

### Domain annotation and rearrangements

We annotated the filtered proteomes with PfamScan v1.6 and Pfam database v33.1:

``` pfam_scan.pl -fasta ../data/proteomes/proteome_loI_noPG.fa -outfile ../data/pfam/proteome_loI_noPG.dom ```

The proteome annotation files from the outgroups (data/proteomes/Elan_pep_noPG_loI.fa, data/proteomes/Focc_pep_noPG_loI.fa) were combined into a single file data/pfam/outgroup.dom

Domain rearrangements were reconstructed with DomRates: https://domainworld.uni-muenster.de/programs/domrates/index.html

``` domRates -t data/tree.nwk -a data/pfam/ -g outgroup -e .dom -o rates.txt -s statistics.txt ```

Only exact solutions were taken in the further analyses.

### Gene Ontology enrichment analysis

* Domain level

To create a GO universe for all the rearranged domains and extract the nodes of interest for the enrichment analysis, the following script was used: scripts/domain2topGo.py

First, we downloaded the pfam2go database

``` wget http://current.geneontology.org/ontology/external2go/pfam2go ```

Then, ran the script domain2topGo.py with the DomRates output and included a node of interest for the enrichment analysis (```-n 7``` in the example):

``` python3 scripts/domain2topGo.py -p data/pfam/ -g pfam2go -s analysis/domrates/statistics_epd.txt -n 7 -c 7node.goi -o Blattodea_universe.go ```

To perform a GO enrichment analysis, we ran the following script (replacing the second argument with genes of interests from different nodes)

1. Universe
2. Genes of interest
3. Algorithm
4. Output

``` scripts/analyseGO_mod.R Blattodea_universe.go 7node.goi weight 7node_weight01 ```

* Gene level

To extract genes with domain rearrangements, we used the following script (specifying a node of interest with the option ``` -n ```):

``` scripts/scripts/extract_DomRates_arrangement_seqs.py -f analysis/domrates/statistics_epd.txt -d data/pfam/ -x dom ```

Then we mapped GO terms from pfam2go and performed GO enrichment analysis like described above. The resulting tables were visualised with the script ``` scripts/GO_heatmaps.R ```

### Differential gene and exon expression analysis

The data for differential gene and exon expression is located in data/ExprData/. A gene was considered to be differentially expressed if adjusted p-value (padj) in at least one comparison is lower than 0.05. The numbers of biased and un-biased genes and rearranged and non-rearranged genes were compared with chi-square test (script ``` scripts/DiffExonExprEnrichment.R ``` for alternatively spliced genes, script ``` scripts/CasteBiasedGenesEnrichment.R ``` for differentially expressed genes).

### Phylostratigraphy analysis

The gene age was calculated with the following program: https://github.com/AlexGa/Phylostratigraphy. The resulting tables used later on are stored in the directory analysis/Phylostrat/

### Rearrangement probability model

The logistic regressions were performed with the following script: ``` scripts/LogitRegrModels.R ```. It considers the following parameters:

1. Gene length (analysis/gene_lengths/)
2. Number of domains per gene (analysis/dom_numbers/)
3. Caste-biased expression (data/ExprData/)
4. Gene age (analysis/Phylostrat/)

The figures illustrating these parameters were generated with the script ``` scripts/figs_and_mods.R ```

### Methylation analysis

First, we extracted the DNA sequences of all rearranged domains with the ``` extractDomains ``` program from DW-Helper (https://domainworld.uni-muenster.de/programs/dw-helper/index.html)

``` extractDomains -d data/pfam/Bger_pep_noPG_loI.dom -i bger_cds.fa -o bger_alldom.fa -D  ```

* Note: The -D option was available only in the development version of ```extractDomains``` at the time when the analysis was performed

The metric CpGoe was used to estimate the methylation level of genes and was calculated with the script ``` scripts/calculate_CpGoe.py  ```. The comparisons of CpGoe between rearranged and non-rearranged genes was performed with the scripts ``` CpGoe_comparison_roaches.R ``` and ``` CpGoe_comparison.R ```.

The methylation data for M. natalensis was taken from the study Harrison et al. 2022 (stored in data/mnat_methylation/) and analysed with the script ``` Methylation_mnat.R  ```

### Transposable elements proximity analysis

To see whether rearranged genes have differences in TE content in the near proximity, we used data from a previous study (Harrison et al. 2018, tables with TE counts per gene are in data/TEproximity/)

Script for the analysis: ``` scripts/TEproximity.R ```

## Contact us

Authors: Alina A. Mikhailova, Elias Dohmen, Mark C. Harrison. 

In case of questions, please contact the first author via email: amikhail[at]uni-muenster.de.


