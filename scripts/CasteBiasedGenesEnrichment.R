library(dplyr)
library(effectsize)
library(pwr)

# Q: is there an over-representation of caste-based genes in the genes with domain rearrangements?

# Csec: total biased and unbiased genes for WF, WM or MF
csec = chisq.test(matrix(c(44,55,4787,13274),nrow = 2))
# data:  matrix(c(44, 55, 4787, 13274), nrow = 2)
# X-squared = 15.323, df = 1, p-value = 9.059e-05

effectsize(csec, type = 'cohens_w')
# Cramer's V (adj.) |       95% CI
# --------------------------------
# 0.03              | [0.02, 1.00]

pwr.chisq.test(w=0.03,
               N = 18160, # total number of observations
               df = 1,
               sig.level=0.05)

# pwr.chisq.test(w=0.03,
#                N = NULL, # total number of observations
#                df = 1,
#                sig.level=0.05,
#                power = 0.98)

# Znev: total biased and unbiased genes for WF, WM or MF
znev = chisq.test(matrix(c(63,42,8657,6695),nrow = 2))
# data:  matrix(c(63, 42, 8657, 6695), nrow = 2)
# X-squared = 0.41565, df = 1, p-value = 0.5191

# Mnat: total biased and unbiased genes for WF, WM or MF
mnat = chisq.test(matrix(c(39,13,6207,9878),nrow = 2))
# data:  matrix(c(39, 13, 6207, 9878), nrow = 2) 
# X-squared = 27.451, df = 1, p-value = 1.611e-07

effectsize(mnat)
# Cramer's V (adj.) |       95% CI
# --------------------------------
# 0.04              | [0.03, 1.00]

pwr.chisq.test(w=0.04,
               N = 16137, # total number of observations
               df = 1,
               sig.level=0.05)

# Rspe: total biased and unbiased genes for W and A (reproductive).
chisq.test(matrix(c(7,66,1476,11250),nrow = 2))
# data:  matrix(c(7, 66, 1476, 11250), nrow = 2)
# X-squared = 0.12353, df = 1, p-value = 0.7252


effectsize(rspe, type="cohens_w")

pwr.chisq.test(w=0.00,
               N = 13948, # total number of observations
               df = 1,
               sig.level=0.05)

