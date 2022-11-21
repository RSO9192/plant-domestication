library(mvabund)
library(lattice)
library(gllvm)
library(corrplot) 
library(phyloseq)
library(tidyverse)
library(metagMisc)
library(metacoder)
library(parallel)
library(readxl)
library(metacoder)
library(spaa)
library(gllvm)

# loading the phyloseq object. 

ps.sv <- readRDS("../../data/Phyloseq_data/ps_sv_PVC.rds")

# Removing SVs that based on 1% prevalence and 20 reads abundance
ps.sv <- phyloseq_filter_prevalence(ps.sv, prev.trh = 0.01, abund.trh = 20, threshold_condition = "OR", abund.type = "total")

# reading in sample data and otu table
samp.data <- data.frame(sample_data(ps.sv))
otu.table <- t(as.data.frame(otu_table(ps.sv)))

# making otu table as mv abund object
otu.MV <- mvabund(otu.table)

#fitting the null model
ft_reveg0=manyglm(otu.MV ~1, data=samp.data)

# offset as difference between total counts peer sample and the NULL model. This account for different sequencing depths
QDoffset = log(rowSums(otu.MV)) - log(rowSums(predict.manyglm(ft_reveg0)))

# Model selection. We do not consider Flowering_time (FT) and seed weight as these variables have not been directly measured 

samp.data2 <- samp.data %>% 
  dplyr::select(Subpopulation, Status, Dom_event,   12:24, -Weight,)

devs = rep(NA,ncol(samp.data2))
names(devs) = colnames(samp.data2)

# loop to calculate the AIC for each variable
for (iVar in 1:ncol(samp.data2)) {
  glmi = manyglm(otu.MV~samp.data2[,iVar],offset = QDoffset,data = samp.data2)
  devs[iVar] = glmi$AICsum
}


# checking the AIC to find the best model

aic <- data.frame(aic=devs, var=names(devs))

aic %>% 
  mutate(var = fct_reorder(var, desc(aic))) %>%
  ggplot(., aes(aic, var))+
  geom_bar(stat = "identity", fill="red", color="black", alpha=0.5)+
  xlab("AIC")+
  ylab("Variable")+
  theme(axis.text.x = element_text(angle=45, vjust=0.5))+
  coord_cartesian(xlim = c(70000, 80000)) 

# the best variable for plant phenotypes is Ca. Now we refit the model to consider Regeneration site and Collection year. 
# We standardize the quantitative variables first. 

# Scaling Ca and Collection year, considering reg site and collection year. 

samp.data3 <- samp.data %>% 
  mutate_at(c("Ca","C.year"), ~(scale(.) %>% as.vector))

glm_best2 = manyglm(otu.MV~Ca+C.year+Regeneration_site,offset = QDoffset,data = samp.data3)

# run from cluster
anova_best2 <- mclapply(1:120, function(x) {anova(glm_best2, nBoot=1,test="LR")}, mc.cores = 120)

anova_PVC_120 <- readRDS("../../data/mvabund_data/anova_vulgaris_120.rds")

K=120 # cores
N=1  #bootstrapping

# p-value for Ca

pi_Ca <- map_dbl(anova_PVC_120, ~ .x$table[2,4])
piMean_Ca = mean(pi_Ca)
(p_Ca = piMean_Ca + (piMean_Ca-1)*(K-1)/(K*N+1)) # Ca is significant. p.value=0.008

# p-value for Sampling time

pi_time <- map_dbl(anova_PVC_120, ~ .x$table[3,4])
piMean_time = mean(pi_time)
(p_time = piMean_time + (piMean_time-1)*(K-1)/(K*N+1)) # Sampling time is not significant. P.value = 0.13

# p-value for Regeneratino site
pi_site <- map_dbl(anova_PVC_120, ~ .x$table[4,4])
piMean_site = mean(pi_site)
(p_site = piMean_site + (piMean_site-1)*(K-1)/(K*N+1)) # Regeneration site is not significant. P.value = 0.13


# Therefore, the minimum adeguate model contains Ca only. 
glm_best3 = manyglm(otu.MV~Ca,offset = QDoffset,data = samp.data)
 
# run this from cluster
anova_PVC_1080 <- mclapply(1:120, function(x) {anova(glm_best2, nBoot=9,test="LR")}, mc.cores = 120)

# indicator species 
out <- anova(glm_best3, nBoot=1,test="LR", p.uni = "adjust") # number of bootstrapping affects p.value only, not deviance explained

# 65% of the model is explaiend by 300 SVs
(sum(sort(out$uni.test[2,], decreasing = TRUE)[1:300]))/out$table[2,3]*100

# extracting coefficients for Ca

coef <- data.frame(Ca=glm_best3$coefficients[2,], SV=names(glm_best3$coefficients[2,]))

# attaching to the tax table the coef for Calcium. Because there were a lot of rare SVS it does not make sense to visualize the 
# coefficients for those SVs as you cannot really fit a model if you have one observation only. For this reason I will visualize the results 
# for the SVs present in at least 5% of the samples.

ps.sv <- readRDS("../../data/Phyloseq_data/ps_sv_PVC.rds")

ps.filt <- phyloseq_filter_prevalence(ps.sv, prev.trh = 0.05, abund.trh = NULL, threshold_condition = "OR", abund.type = "total")

# taxas filtered at 5% prevalence
tax.table <- data.frame(tax_table(ps.filt))

tax.table$SV <- rownames(tax.table)

# combining coefficients for taxa present at 5% prevalence
taxmtab.coef <- left_join(tax.table, coef)

# visualizing relative importance of positive and negative slopes for filtered data at 5%

taxmtab.coef %>% 
  filter(Phylum%in%c("Proteobacteria", "Acidobacteriota", "Actinobacteriota", "Bacteroidota", "Deinococcota", "Firmicutes")) %>% 
  mutate(Ca=ifelse(Ca<0,"Negative", "Positive")) %>% 
  group_by(Phylum, Ca) %>% 
  summarise(amount=n()) %>% 
  ungroup() %>% 
  group_by(Phylum) %>% 
  mutate(tot=sum(amount), r.a=(amount/tot)*100) %>% 
  ggplot()+
  geom_bar(aes(Phylum,r.a,fill=Ca),stat="identity", color="black")+
  ylab("Relative importance of coefficients")+
  theme(axis.text.x = element_text(
    size = 12, angle = 45, hjust = 1), axis.text.y=element_blank(),
    axis.ticks.y=element_blank())
