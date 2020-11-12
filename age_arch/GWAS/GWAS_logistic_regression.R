f = '/Users/fongsl/Desktop/FANTOM_x_gwas19_ldex.bed'
install.packages('aod')
library(aod)
install.packages('ggplot2')
library(ggplot)

df <- read.csv(f, sep= '\t', header = FALSE)
keep <-df[ , c("V1", "V2", "V3","V4",  "V5",  "V6",  "V7", "V8", "V16")] 




colnames(keep) <-  c("chr_enh", "start_enh", "end_enh", "old_len", "core_remodeling",
                  "mrca_2", "datatype", "count_overlap", "gwas_overlap")

summary(keep)

sapply(keep, sd)
head(keep)

keep$core_remodeling<- as.factor(keep$core_remodeling)
keep$gwas_overlap <- as.factor(keep$gwas_overlap)
keep$mrca_2<- as.factor(keep$mrca_2)

rawlogit <- glm(gwas_overlap ~ count_overlap + mrca_2 + core_remodeling + old_len, data = keep, family = binomial(("logit")))
summary(rawlogit)

keep$wts  <- ifelse(keep$gwas_overlap ==1, 1, (1281/30218))
head(keep)


wlogit <- glm(gwas_overlap ~ count_overlap + mrca_2 + core_remodeling + old_len, data = keep, weights = wts, family = binomial(("logit")))
summary(wlogit)
wlogit_interaction <- glm(gwas_overlap ~ count_overlap + mrca_2 + core_remodeling + old_len + count_overlap*core_remodeling, data = keep, weights = wts, family = binomial(("logit")))
summary(wlogit_interaction)
