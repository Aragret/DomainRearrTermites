library(ggplot2)


# comparing dogma scores between different groups (cockroaches vs termites, species with and without true workers)


dogma = c(75.72, 86.16, 59.70, 68.82, 84.61, 90.01, 83.64, 93.88, 84.97, 75.66)
species = c("elan", 'focc', 'bger', 'dpun', 'pame', 'znev', 'csec', 'rspe', 'cfor', 'mnat')
termites = factor(c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1))
workers = factor(c(0, 0, 0, 0, 0, 0, 0, 1, 1, 1))

df = data.frame(sp = species, dogma = dogma, termites = termites, workers = workers)

summary(aov(dogma ~ termites, df))
# Df Sum Sq Mean Sq F value Pr(>F)
# termites     1  282.5  282.49    3.32  0.106
# Residuals    8  680.7   85.09

summary(aov(dogma ~ workers, df))
# Df Sum Sq Mean Sq F value Pr(>F)
# workers      1   87.5   87.55     0.8  0.397
# Residuals    8  875.7  109.46

summary(aov(dogma ~ termites + workers, df))
# Df Sum Sq Mean Sq F value Pr(>F)
# termites     1  282.5  282.49   2.925  0.131
# workers      1    4.7    4.74   0.049  0.831
# Residuals    7  676.0   96.57 

ggplot(df, aes(termites, dogma)) + 
  geom_boxplot()

ggplot(df, aes(workers, dogma)) + 
  geom_boxplot()
