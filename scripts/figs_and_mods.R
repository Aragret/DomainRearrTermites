install.packages("ggsci")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggstatsplot)
library(ggsci)


# install.packages("ggiraphExtra")
library(ggiraphExtra)

install.packages("predict3d")
library()

### Input

LengthsTable = "Mnat_gene_lengths"
DiffExprData = "bigtable_mnat_newest"
OutputFig = "Mnat_lengths.pdf"
# RA = "analysis/domrates/termites_arrangements_exact_genes"
RA = "analysis/domrates/advanced_sociality_arrangements_exact_genes"
Domains = "analysis/dom_numbers/Mnat_dom_numbers"
AS = "data/ExprData/all_as_01"

lengths = read.table(LengthsTable, sep = "\t", header = TRUE)
rna = read.table(DiffExprData, sep = "\t", header = TRUE)
RA = read.table(RA)
Doms = read.table(Domains, sep = ' ')
AS = read.table(AS)

names(lengths)[1] = "gene_id"
names(RA) = "gene_id"
RA$RA = 1
names(Doms) = c("Doms", "gene_id")
names(AS) = "gene_id"
AS$AS = 1

# for csec
# lengths = separate(data = lengths, GeneID, c('gene_id', NA), '-')

df = right_join(rna, lengths) %>%
  left_join(RA) %>%
  left_join(AS) %>%
  mutate(
    Biased = as.factor(case_when(
      padj_WF < 0.05 | padj_WM < 0.05 | padj_MF < 0.05 ~ 1,
      padj_WF >= 0.05 | padj_WM >= 0.05 | padj_MF >= 0.05 ~ 0
    )),
    RA = as.factor(case_when(is.na(RA) ~ 0,
                   RA == 1 ~ 1)),
    AS = as.factor(case_when(is.na(AS) ~ 0,
                   AS == 1 ~ 1)),
  ) %>%
  left_join(mnat_bigdf[, c("gene_id", "PS1")]) %>%
  left_join(Doms)

# df$Doms[is.na(df$Doms)] = 0

pdf(OutputFig)

print(ggbetweenstats(
  data = df,
  x = Biased,
  y = Length,
  type = "nonparametric",
  title = "Distribution of gene lengths across caste-biased and unbiased genes"
))

print(ggbetweenstats(
  data = df,
  x = AS,
  y = Length,
  type = "nonparametric",
  title = "Distribution of gene lengths across genes with caste-biased exon expression"
))

print(ggbetweenstats(
  data = df,
  x = RA,
  y = Length,
  type = "nonparametric",
  title = "Distribution of gene lengths across rearranged and non-rearranged genes"
))

print(ggbetweenstats(
  data = df,
  x = Biased,
  y = Doms,
  type = "nonparametric",
  title = "Distribution of domain numbers across caste-biased and unbiased genes"
))

print(ggbetweenstats(
  data = df,
  x = AS,
  y = Doms,
  type = "nonparametric",
  title = "Distribution of domain numbers across genes with caste-biased exon expression"
))

print(ggbetweenstats(
  data = df,
  x = RA,
  y = Doms,
  type = "nonparametric",
  title = "Distribution of domain numbers across rearranged and non-rearranged genes"
))

dev.off()

# ggplot(df[!is.na(df$Biased),], aes(RA, Doms, fill = Biased)) +
#   geom_violin()

summary(glm(RA ~ Doms + Biased + Doms*Biased, data = df, family = "binomial"))

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -4.81678    0.16155 -29.816  < 2e-16 ***
#   Doms          0.10556    0.04411   2.393  0.01670 *  
#   Biased1       0.52184    0.19855   2.628  0.00858 ** 
#   Doms:Biased1 -0.06580    0.04950  -1.329  0.18379    

# without genes with 0 domains

# (Intercept)  -4.12418    0.21343 -19.324   <2e-16 ***
#   Doms         -0.03818    0.08920  -0.428    0.669    
# Biased1       0.21360    0.25320   0.844    0.399    
# Doms:Biased1  0.03154    0.09832   0.321    0.748    

summary(glm(RA ~ Doms + AS + Doms*AS, data = df, family = "binomial"))

# (Intercept) -4.224789   0.128216 -32.951   <2e-16 ***
#   Doms        -0.003184   0.050364  -0.063   0.9496    
# AS1          0.615614   0.253585   2.428   0.0152 *  
#   Doms:AS1     0.001374   0.068832   0.020   0.9841  

summary(glm(RA ~ Doms + AS + Length + Doms*AS + Length*AS, data = df, family = "binomial"))

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.3491751  0.1399142 -31.085  < 2e-16 ***
#   Doms        -0.0881769  0.0652337  -1.352  0.17647    
# AS1          0.6865826  0.2671693   2.570  0.01017 *  
#   Length       0.0006049  0.0001910   3.167  0.00154 ** 
#   Doms:AS1     0.0617955  0.0886990   0.697  0.48600    
# AS1:Length  -0.0004692  0.0002632  -1.783  0.07464 .  

summary(glm(RA ~ Doms + Biased + Length, data = df, family = "binomial"))

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.2168787  0.1590674 -26.510  < 2e-16 ***
#   Doms        -0.0778291  0.0492539  -1.580  0.11407    
# Biased1      0.2418020  0.1785971   1.354  0.17577    
# Length       0.0003455  0.0001184   2.918  0.00352 ** 

summary(glm(RA ~ Doms + Biased + Length + Biased*Doms + Biased*Length, data = df, family = "binomial"))

# Estimate Std. Error z value Pr(>|z|)    
# (Intercept)    -4.189e+00  2.233e-01 -18.759   <2e-16 ***
#   Doms           -8.373e-02  1.055e-01  -0.794    0.427    
# Biased1         2.083e-01  2.624e-01   0.794    0.427    
# Length          3.120e-04  2.958e-04   1.055    0.292    
# Doms:Biased1    7.349e-03  1.190e-01   0.062    0.951    
# Biased1:Length  3.809e-05  3.232e-04   0.118    0.906  


### compare different models 

# RA 0 - 7871, RA 1 - 140

m0 = glm(RA ~ 1, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")

Biased = glm(RA ~ Biased, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")
summary(Biased)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept)  -4.1942     0.1425 -29.436   <2e-16 ***
#   Biased1       0.2698     0.1779   1.517    0.129 

anova(Biased, m0, test="LR")
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1      8009     1408.3                     
# 2      8010     1410.7 -1   -2.353    0.125

AS = glm(RA ~ AS, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")
summary(AS)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.13120    0.09745 -42.394  < 2e-16 ***
#   AS1          0.52933    0.20156   2.626  0.00863 ** 

anova(AS, m0, test="LR")
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
# 1      8009     1404.4                       
# 2      8010     1410.7 -1  -6.2781  0.01222 *

Length = glm(RA ~ Length, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")
summary(Length)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.142e+00  1.005e-01 -41.197   <2e-16 ***
#   Length       2.016e-04  8.614e-05   2.341   0.0193 * 

anova(Length, m0, test="LR")
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
# 1      8009     1406.9                       
# 2      8010     1410.7 -1  -3.8212  0.05061 .

Doms = glm(RA ~ Doms, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")
summary(Doms)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.00798    0.11409 -35.130   <2e-16 ***
#   Doms        -0.01074    0.03872  -0.277    0.781   

anova(Doms, m0, test="LR")
# Resid. Df Resid. Dev Df  Deviance Pr(>Chi)
# 1      8009     1410.6                      
# 2      8010     1410.7 -1 -0.082294   0.7742

Biased_Length = glm(RA ~ Biased + Length, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")
summary(Biased_Length)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.285e+00  1.491e-01 -28.743   <2e-16 ***
#   Biased1      2.440e-01  1.786e-01   1.366   0.1718    
# Length       1.904e-04  8.739e-05   2.179   0.0293 *  

anova(Biased_Length, m0, test="LR")
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
# 1      8008     1405.0                       
# 2      8010     1410.7 -2  -5.7251  0.05712 .

AS_Length = glm(RA ~ AS + Length, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")
summary(AS_Length)
# Coefficients:
#   Estimate Std. Error z value Pr(>|z|)    
# (Intercept) -4.202e+00  1.076e-01 -39.042   <2e-16 ***
#   AS1          4.564e-01  2.093e-01   2.180   0.0293 *  
#   Length       1.506e-04  9.438e-05   1.596   0.1105 

anova(AS_Length, m0, test="LR")
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)  
# 1      8008     1402.5                       
# 2      8010     1410.7 -2  -8.2257  0.01636 *

anova(AS_Length, AS, test="LR")
# Resid. Df Resid. Dev Df Deviance Pr(>Chi)
# 1      8008     1402.5                     
# 2      8009     1404.4 -1  -1.9477   0.1628

m1 = glm(RA ~ Biased + AS, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")

m2 = glm(RA ~ Biased + AS + Doms + Biased*Doms + AS*Doms, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")

m3 = glm(RA ~ Biased + AS + Length, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")

m4 = glm(RA ~ Biased + AS + Length + Biased*Length + AS*Length, data = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),], family = "binomial")


anova(m2, m1, test="LR") # not better

anova(m3, m1, test="LR") # not better

anova(m4, m1, test="LR") # not better

### with gene age

summary(glm(RA ~ AS + Biased + Length + PS1, data = df, family = "binomial"))
# AS1          4.829e-01  2.108e-01   2.291    0.022 *  
#   Biased1      2.433e-01  1.793e-01   1.357    0.175    
# Length       1.458e-04  9.453e-05   1.542    0.123    
# PS11         1.701e+00  3.330e-01   5.108 3.25e-07 ***

summary(glm(RA ~ AS + Biased + Length + PS1 + Doms, data = df, family = "binomial"))
# AS1          0.4064880  0.2106156   1.930   0.0536 .  
# Biased1      0.2081467  0.1789506   1.163   0.2448    
# Length       0.0002898  0.0001254   2.310   0.0209 *  
#   PS11         0.5903882  0.3331358   1.772   0.0764 .  
# Doms        -0.0869276  0.0522593  -1.663   0.0962 .  

summary(glm(RA ~ AS * Biased * Length * PS1 * Doms, data = df, family = "binomial"))

ggplot(df, aes(x=Length, y=as.numeric(as.character(RA)), col=AS)) + 
  geom_point(alpha=0.3) +
  stat_smooth(method="glm", color="green", se=FALSE, 
              method.args = list(family=binomial))

ggplot(df, aes(x=Length, y=as.numeric(as.character(RA)))) + 
  geom_point(alpha=0.3) +
  stat_smooth(method="glm", se=FALSE, 
              method.args = list(family=binomial)) +
  facet_grid(~ AS)

ggplot(df, aes(x=Length, y=as.numeric(as.character(RA)))) + 
  geom_point(alpha=0.3) +
  stat_smooth(method="glm", color="green", se=FALSE, 
              method.args = list(family=binomial)) +
  facet_grid(~ Biased)

test = df[!is.na(df$Biased) & !is.na(df$AS) & !is.na(df$Doms),]
test$RA = as.numeric(as.character(test$RA))
test_m = glm(RA ~ Length + AS, test, family = "binomial")

test_m2 = glm(RA ~ Length, test, family = "binomial")

ggPredict(test_m)
ggPredict(test_m2)


### csec
table(df[df$RA == 1,]$Biased)
# 0  1 
# 53 44 
table(df[df$RA == 0,]$Biased)
# 0    1 
# 9097 4787 

44/(44+53) # 0.4536082
4787/(4787+9097) # 0.3447854

chisq.test(matrix(c(53,44,9097,4787),nrow = 2))
# X-squared = 4.5746, df = 1, p-value = 0.03245

### znev
table(df[df$RA == 1,]$Biased)
# 0  1 
# 40 63 
table(df[df$RA == 0,]$Biased)
# 0    1 
# 5271 8657 

63/(63+40) # 0.6116505
8657/(8657+5271) # 0.6215537

chisq.test(matrix(c(40,63,5271,8657),nrow = 2))
# X-squared = 0.010922, df = 1, p-value = 0.9168

### mnat origin of termites

table(df[df$RA == 1,]$Biased)
# 0  1 
# 13 39

table(df[df$RA == 0,]$Biased)
# 0    1 
# 5360 6207

39/(39+13) # 0.75

6207/(5360+6207) # 0.5366128

chisq.test(matrix(c(13,39,5360,6207),nrow = 2))
# X-squared = 8.6433, df = 1, p-value = 0.003283

### mnat origin of true workers

table(df[df$RA == 1,]$Biased)
# 0  1 
# 50 90 

table(df[df$RA == 0,]$Biased)
# 0    1 
# 5323 6156

90/(90+50) # 0.6428571
6156/(6156+5323) # 0.5362836

chisq.test(matrix(c(50,90,5323,6156),nrow = 2))
# X-squared = 5.898, df = 1, p-value = 0.01516

