##-----------Installing the necessary packages if not installed
if (!require("readr")) install.packages("plotly")
if (!require("DescTools")) install.packages("brunnermunzel")
if (!require("ggplot2")) install.packages("readr")

##-----------Loading the libraries required in this R code
library(readr)
library(DescTools)
library(ggplot2)

##-----------Importing the 'Total occurances per disease.csv'file as df from the 
##-----------python code output
df <- read_csv("Total occurances per disease.csv")
df <- df[df$Status=='reviewed',]
df$geneage <- factor(df$`Gene age`, labels = c("recent", "young", "old", "old", "old", "archaic", "archaic", "primordial"),
                                    levels = c("Mammalia", "Vertebrata", "Eumetazoa", "Opisthokonta", "Eukaryota", "Euk_Archaea", "Euk+Bac", "Cellular_organisms"), ordered = T)
df$mutual <- cut(df$`Mutual information estimation (in pits)`, c(700, 3000, Inf), right = T, ordered_result = T)
df$disease_or_nondisease <- factor(!is.na(df$`Involvement in disease (source: UniProt)`), ordered = T)


# disease or no disease tests

a <- chisq.test(df$disease_or_nondisease, df$`Number of interactions`)
b <- CramerV(df$disease_or_nondisease, df$`Number of interactions`)
cat('test statistics for correlation between diseases boolean and number of interaction
    are:\nchi-square:\t', a$statistic,"\np-value:\t", a$p.value, "\nCrammers's V:\t", b)

c <- chisq.test(df$disease_or_nondisease, df$`Protein length`)
d <- CramerV(df$disease_or_nondisease, df$`Protein length`)
cat('test statistics for correlation between diseases boolean and Protein length
    are:\nchi-square:\t', c$statistic,"\np-value:\t", c$p.value, "\nCrammers's V:\t", d)

e1 <- chisq.test(df$disease_or_nondisease, df$mutual)
e2 <- CramerV(df$disease_or_nondisease, df$mutual)
cat('test statistics for correlation between diseases boolean and mutual information ranks
    are:\nchi-square:\t', e1$statistic,"\np-value:\t", e1$p.value, "\nCrammers's V:\t", e2)

df1 <- df[!is.na(df$geneage),]
f <- KendallTauB(df1$disease_or_nondisease, df1$geneage)
cat('test statistics for correlation between diseases boolean and gene ages
    are:\nKendall TauB:\t', f)
    
# disease cateogory tests

df$disease_occurance_cateogory <- cut(df$`Total occurance value`, c(0, 0.000001, 0.0001, 1), right = T, ordered_result = T)
df2 <- df[!is.na(df$disease_occurance_cateogory),]

g <- chisq.test(df2$disease_occurance_cateogory, df2$`Number of interactions`)
h <- CramerV(df2$disease_occurance_cateogory, df2$`Number of interactions`)
cat('test statistics for correlation between disease cateogories and number of interaction
    are:\nchi-square:\t', g$statistic,"\np-value:\t", g$p.value, "\nCrammers's V:\t", h)

i <- chisq.test(df2$disease_occurance_cateogory, df2$`Protein length`)
j <- CramerV(df2$disease_occurance_cateogory, df2$`Protein length`)
cat('test statistics for correlation between disease cateogories and Protein length
    are:\nchi-square:\t', i$statistic,"\np-value:\t", i$p.value, "\nCrammers's V:\t", j)

k1 <- chisq.test(df2$disease_occurance_cateogory, df2$mutual)
k2 <- CramerV(df2$disease_occurance_cateogory, df2$mutual)
cat('test statistics for correlation between disease cateogories and mutual information ranks
    are:\nchi-square:\t', k1$statistic,"\np-value:\t", k1$p.value, "\nCrammers's V:\t", k2)

df3 <- df2[!is.na(df2$geneage),]
l <- KendallTauB(df3$disease_or_nondisease, df3$geneage)
cat('test statistics for correlation between diseases boolean and gene ages
    are:\nKendall TauB:\t', l)



ggplot(data = df1, mapping = aes(geneage, col = disease_or_nondisease, fill=disease_or_nondisease)) +
  geom_density(adjust=.211, alpha=.15)

ggplot(df1, aes(geneage, col=disease_or_nondisease)) +
  geom_freqpoly(binwidth = 40) +
  xlim(1, 200)

#results

cat('test statistics for correlation between diseases boolean and number of interaction
    are:\nchi-square:\t', a$statistic,"\np-value:\t", a$p.value, "\nCrammers's V:\t", b)


cat('test statistics for correlation between diseases boolean and Protein length
    are:\nchi-square:\t', c$statistic,"\np-value:\t", c$p.value, "\nCrammers's V:\t", d)

cat('test statistics for correlation between diseases boolean and mutual information ranks
    are:\nchi-square:\t', e1$statistic,"\np-value:\t", e1$p.value, "\nCrammers's V:\t", e2)

cat('test statistics for correlation between diseases boolean and gene ages
    are:\nKendall TauB:\t', f)

cat('test statistics for correlation between disease cateogories and number of interaction
    are:\nchi-square:\t', g$statistic,"\np-value:\t", g$p.value, "\nCrammers's V:\t", h)

cat('test statistics for correlation between disease cateogories and Protein length
    are:\nchi-square:\t', i$statistic,"\np-value:\t", i$p.value, "\nCrammers's V:\t", j)

cat('test statistics for correlation between disease cateogories and mutual information ranks
    are:\nchi-square:\t', k1$statistic,"\np-value:\t", k1$p.value, "\nCrammers's V:\t", k2)

cat('test statistics for correlation between diseases boolean and gene ages
    are:\nKendall TauB:\t', l)


