library(dplyr)
library(tidyr)
library(ggplot2)
library(ggstatsplot)
library(data.table)

### Input

LengthsTable = "analysis/gene_lengths/Rspe_gene_lengths"
DiffExprData = "data/ExprData/bigtable_rspe"
OutputFig = "Rspe_lengths.pdf"
RA = "analysis/domrates/termites_arrangements_exact_genes"
# RA = "analysis/domrates/advanced_sociality_arrangements_exact_genes"
Domains = "analysis/dom_numbers/Rspe_dom_numbers"
PS = read.table("analysis/Phylostrat/Rspe_Age.tsv",
                header = TRUE, sep = "\t")

### Input mnat 

LengthsTable = "analysis/gene_lengths/Mnat_gene_lengths"
DiffExprData = "data/ExprData/bigtable_mnat_newest"
# OutputFig = "Mnat_lengths.pdf"
RA = "analysis/domrates/termites_arrangements_exact_genes"
# RA = "analysis/domrates/advanced_sociality_arrangements_exact_genes"
Domains = "analysis/dom_numbers/domain_stats/Mnat_dom_numbers"
AS = "data/ExprData/all_as_01"
PS = read.table("analysis/Phylostrat/Mnat_Age.tsv",
                header = TRUE, sep = "\t")

OutputFig_geneAge = "GeneAge_mnat_workers.pdf"

lengths = read.table(LengthsTable, sep = "\t", header = TRUE)
rna = read.table(DiffExprData, sep = "\t", header = TRUE)

#for rspe
# rna = read.table(DiffExprData, sep = ",", header = TRUE)
# names(rna)[1] = "gene_id"

RA = read.table(RA)
Doms = read.table(Domains, sep = ' ')

# mnat
# AS = read.table(AS)
# names(AS) = "gene_id"
# AS$AS = 1

names(lengths)[1] = "gene_id"
names(RA) = "gene_id"
RA$RA = 1
names(Doms) = c("Doms", "gene_id")
names(PS)[2] = "gene_id"



# csec
# lengths = separate(data = lengths, gene_id, c('gene_id', NA), '-')
# PS = separate(data = PS, gene_id, c('gene_id', NA), '-')
# Doms = separate(data = Doms, gene_id, c('gene_id', NA), '-')

# rspe
# lengths = separate(data = lengths, gene_id, c('gene_id', NA), '-')
# Doms = separate(data = Doms, gene_id, c('gene_id', NA), '-')
# RA$gene_id = gsub("-PA", "", RA$gene_id)
# PS$gene_id = gsub("-PA", "", PS$gene_id)



df = right_join(rna, lengths) %>% # -6 for mnat, -12 for others
  left_join(RA) %>%
  mutate(
    Biased = as.factor(case_when(
      padj < 0.05 ~ 1,
      padj >= 0.05 ~ 0)),
    # Biased = as.factor(case_when(
    #   padj_WF < 0.05 | padj_WM < 0.05 | padj_MF < 0.05 ~ 1,
    #   padj_WF >= 0.05 | padj_WM >= 0.05 | padj_MF >= 0.05 ~ 0)),
    RA = as.factor(case_when(is.na(RA) ~ 0,
                             RA == 1 ~ 1))) %>%
  left_join(Doms) %>% left_join(PS) %>%
  mutate(PS_groups = as.factor(case_when(
    PS == 1 ~ "PS1_Cellular_Organisms",
    PS >= 2 & PS <= 9 ~ "PS2_Eukaryota",
    PS >= 10 & PS <= 16 ~ "PS3_Arthropoda",
    PS >= 17 & PS <= 18 ~ "PS4_Neoptera",
    PS == 19 ~ "PS5_Dictyoptera",
    PS >= 20 & PS <= 22 ~ "PS6_Termitoidea",
    PS >= 23 ~ "PS7_Species_specific"
  ))) %>%
  # mnat
  left_join(AS) %>% 
  mutate(
    AS = as.factor(case_when(is.na(AS) ~ 0,
                             AS == 1 ~ 1)),
  )

df_models = df[!is.na(df$Biased) & !is.na(df$Doms) & !is.na(df$PS_groups),]


# make PS dummy

# df_PS_groups = left_join(dcast(df_models, gene_id ~ PS_groups,
#       fun.aggregate = function(x) 1L,
#       fill = 0L), df_models[, c("gene_id", "RA", "Length", "Doms", "Biased")])

df_PS_groups = df_models %>%
  mutate(
    PS1_Cellular_Organisms = as.factor(case_when(
      PS_groups == "PS1_Cellular_Organisms" ~ 1,
      PS_groups != "PS1_Cellular_Organisms" ~ 0,
    )),
    PS2_Eukaryota = as.factor(case_when(
      PS_groups %in% c("PS1_Cellular_Organisms") ~ 0,
      !(PS_groups %in% c("PS1_Cellular_Organisms")) ~ 1
    )),
    PS3_Arthropoda = as.factor(case_when(
      PS_groups %in% c("PS1_Cellular_Organisms", "PS2_Eukaryota") ~ 0,
      !(PS_groups %in% c("PS1_Cellular_Organisms", "PS2_Eukaryota")) ~ 1
    )),
    PS4_Neoptera = as.factor(case_when(
      PS_groups %in% c("PS1_Cellular_Organisms", "PS2_Eukaryota", "PS3_Arthropoda") ~ 0,
      !(PS_groups %in% c("PS1_Cellular_Organisms", "PS2_Eukaryota", "PS3_Arthropoda")) ~ 1
    )),
    PS5_Dictyoptera = as.factor(case_when(
      PS_groups %in% c("PS1_Cellular_Organisms", "PS2_Eukaryota", "PS3_Arthropoda", "PS4_Neoptera") ~ 0,
      !(PS_groups %in% c("PS1_Cellular_Organisms", "PS2_Eukaryota", "PS3_Arthropoda", "PS4_Neoptera")) ~ 1
    )),
    PS6_Termitoidea = as.factor(case_when(
      PS_groups %in% c("PS1_Cellular_Organisms", "PS2_Eukaryota", "PS3_Arthropoda", "PS4_Neoptera",
                       "PS5_Dictyoptera") ~ 0,
      !(PS_groups %in% c("PS1_Cellular_Organisms", "PS2_Eukaryota", "PS3_Arthropoda", "PS4_Neoptera",
                         "PS5_Dictyoptera")) ~ 1
    )),
    PS7_Species_specific = as.factor(case_when(
      PS_groups %in% c("PS7_Species_specific") ~ 1,
      !(PS_groups %in% c("PS7_Species_specific")) ~ 0
    ))
  )

table(df_PS_groups$PS_groups)

summary(glm(RA ~ PS1_Cellular_Organisms, df_PS_groups, family = "binomial"))

summary(glm(RA ~ PS2_Eukaryota, df_PS_groups, family = "binomial"))

summary(glm(RA ~ PS3_Arthropoda, df_PS_groups, family = "binomial"))

summary(glm(RA ~ PS4_Neoptera, df_PS_groups, family = "binomial"))

summary(glm(RA ~ PS5_Dictyoptera, df_PS_groups, family = "binomial"))

summary(glm(RA ~ PS6_Termitoidea, df_PS_groups, family = "binomial"))

summary(glm(RA ~ PS7_Species_specific, df_PS_groups, family = "binomial"))

summary(glm(RA ~ PS_groups, df_models, family = "binomial"))

table(df_models$RA)

# rspe
# 0    1 
# 9174   71 

# znev 
# 0    1 
# 9417  101 

# csec
# 0     1 
# 10225    96 

# mnat (termites)
# 0    1 
# 7933   52 

pdf(OutputFig_geneAge)

ggplot(df[!is.na(df$PS_groups),], aes(RA, fill=PS_groups)) +
  geom_bar(position = "fill") + 
  theme_classic() +
  scale_fill_simpsons(alpha = 0.9) +
  labs(title = "M. natalensis", xlab = "Rearranged") 

dev.off()


summary(glm(RA ~ Biased + Length + Doms + PS_groups, df_models, family = "binomial"))

summary(glm(RA ~ Biased + Length + PS_groups, df_models, family = "binomial"))

# summary(glm(RA ~ Biased + Length + Doms + PS1_Cellular_Organisms, df_PS_groups, family = "binomial"))

# summary(glm(RA ~ Biased + Length + PS1_Cellular_Organisms, df_PS_groups, family = "binomial"))


summary(glm(RA ~ Biased * Length * PS1_Cellular_Organisms * Doms, df_PS_groups, family = "binomial"))

summary(glm(RA ~ Biased * Length * PS1_Cellular_Organisms, df_PS_groups, family = "binomial"))



summary(glm(RA ~ Biased * Length * PS2_Eukaryota, df_PS_groups, family = "binomial"))

summary(glm(RA ~ Biased * Length * PS3_Arthropoda, df_PS_groups, family = "binomial"))

summary(glm(RA ~ Biased * Length * PS4_Neoptera, df_PS_groups, family = "binomial"))

summary(glm(RA ~ Biased * Length * PS5_Dictyoptera, df_PS_groups, family = "binomial"))

summary(glm(RA ~ Biased * Length * PS6_Termitoidea, df_PS_groups, family = "binomial"))

summary(glm(RA ~ Biased * Length * PS7_Species_specific, df_PS_groups, family = "binomial"))


summary(glm(RA ~ Biased * Length * PS_groups, df_PS_groups, family = "binomial"))

### AS


summary(glm(RA ~ Biased + Length + Doms + PS_groups + AS, df_models, family = "binomial"))

summary(glm(RA ~ Biased + Length + PS_groups + AS, df_models, family = "binomial"))

# summary(glm(RA ~ Biased + Length + Doms + PS1_Cellular_Organisms, df_PS_groups, family = "binomial"))

# summary(glm(RA ~ Biased + Length + PS1_Cellular_Organisms, df_PS_groups, family = "binomial"))

summary(glm(RA ~ Biased * Length * PS1_Cellular_Organisms * AS, df_PS_groups, family = "binomial"))

