# This script was used for the model-based multivariate analysis. Here we show how to perform this analysis for the 
# experiment with P. vulgaris from CGIAR-CIAT (2nd experiment). The same logic applies for the other experiments. PVC in the file 
# name stands for P. vulgaris CIAT

library(mvabund)
library(lattice)
library(corrplot) 
library(phyloseq)
library(tidyverse)
library(metagMisc)
library(parallel)
library(ecoCopula)

# setting working directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# loading the phyloseq object. All phyloseq objects can be downloaded from GitHub.

ps.sv <- readRDS("../data/Phyloseq_data/ps_sv_PVC.rds")

# Removing SVs that based on 1% prevalence and 20 reads abundance. This was applied to all phyloseq objects prior to any statistical analysis.
ps.sv <- phyloseq_filter_prevalence(ps.sv, prev.trh = 0.01, abund.trh = 20, threshold_condition = "OR", abund.type = "total")


# mvabund -----------------------------------------------------------------

# reading in sample data and otu table
samp.data <- data.frame(sample_data(ps.sv))
sv.table <- t(as.data.frame(otu_table(ps.sv)))

# making otu table as mv abund object
sv.MV <- mvabund(sv.table)

#fitting the null model. Default in manyglm is negative bnomial
ft_null=manyglm(sv.MV ~1, data=samp.data)

# offset as difference between total counts per sample and the NULL model. This account for different sequencing depths.
# for more info see (https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13703?af=R)
QDoffset = log(rowSums(sv.MV)) - log(rowSums(predict.manyglm(ft_null)))

# Model selection. We do not consider Flowering_time (FT) and seed weight as these variables have not been directly measured 

samp.data2 <- samp.data %>% 
  dplyr::select(Subpopulation, Status, Dom_event,   12:24, -Weight,) %>% 
  mutate(Dom_event=case_when(Subpopulation=="WA1" ~ "AW",
                             Subpopulation=="M3" ~ "MW",
                             Subpopulation=="M4" ~ "MW",
                             Subpopulation=="WA2" ~ "AW",
                             Subpopulation=="DA1" ~ "AD",
                             Subpopulation=="M2" ~ "MD"))

aic = rep(NA,ncol(samp.data2))
names(aic) = colnames(samp.data2)

# calculating AIC for each variable. COnsidering a quadratic term only when it leads to a lower AIC
for (iVar in 1:ncol(samp.data2)) {
  
  if(is.numeric(samp.data2[,iVar])) {
    
    glms = manyglm(sv.MV~ samp.data2[,iVar], data = samp.data2, offset = QDoffset)
    glms_poly = manyglm(sv.MV~ poly(samp.data2[,iVar],2), data = samp.data2, offset = QDoffset)
    
    if(glms_poly$AICsum <  glms$AICsum) { # inducing quadratic term only when better
      
      aic[iVar] = glms_poly$AICsum 
      
    } else {
      
      aic[iVar] = glms$AICsum 
      
    }
    
    
  } else {
    
    glms = manyglm(sv.MV~ samp.data2[,iVar], data = samp.data2, offset = QDoffset)
    aic[iVar] = glms$AICsum
  }
  
}


# checking the AIC to find the best model

aic <- data.frame(aic=aic, var=names(aic))

aic %>% 
  mutate(var = fct_reorder(var, desc(aic))) %>%
  ggplot(., aes(aic, var))+
  geom_bar(stat = "identity", fill="red", color="black", alpha=0.5)+
  theme_bw(base_size = 11)+
  ylab("AIC")+
  theme(axis.text.x = element_text(angle=45, vjust=0.5))+
  coord_cartesian(xlim = c(72500, 78900)) 

# the best explanatory variable is Ca according to AIC. Now we refit the model to consider Regeneration site and Collection year. 
# We standardize the quantitative variables first. This step is done only for P. vulgaris and P. lunatus grown at CGIAR-CIAT

# Scaling Ca and Collection year, considering reg site and collection year. 

samp.data3 <- samp.data %>% 
  mutate_at(c("Ca","C.year"), ~(scale(.) %>% as.vector))

glm_tot = manyglm(sv.MV~Ca+C.year+Regeneration_site,offset = QDoffset,data = samp.data3)

# run from cluster
anova_glm_tot <- mclapply(1:120, function(x) {anova(glm_tot, nBoot=1,test="LR")}, mc.cores = 120)

# the outout of this analysis can also be downloaded from GitHub
anova_glm_tot <- readRDS("../data/mvabund_data/anova_vulgaris_120.rds")

K=120 # cores
N=1  #bootstrapping

# p-value for Ca

pi_Ca <- map_dbl(anova_glm_tot, ~ .x$table[2,4])
piMean_Ca = mean(pi_Ca)
(p_Ca = piMean_Ca + (piMean_Ca-1)*(K-1)/(K*N+1)) # Ca is significant. p.value=0.008

# p-value for Sampling time

pi_time <- map_dbl(anova_glm_tot , ~ .x$table[3,4])
piMean_time = mean(pi_time)
(p_time = piMean_time + (piMean_time-1)*(K-1)/(K*N+1)) # Sampling time is not significant. P.value = 0.13

# p-value for Regeneratino site
pi_site <- map_dbl(anova_glm_tot , ~ .x$table[4,4])
piMean_site = mean(pi_site)
(p_site = piMean_site + (piMean_site-1)*(K-1)/(K*N+1)) # Regeneration site is not significant. P.value = 0.26


# Therefore, the minimum adequate model contains Ca only. 
glm_best = manyglm(sv.MV~Ca,offset = QDoffset,data = samp.data)
# cannot see clear problems with residuals, so it is ok
plot(glm_best)
# run this from cluster for faster computation
anova_best <- mclapply(1:120, function(x) {anova(glm_best, nBoot=9,test="LR")}, mc.cores = 120)

# extracting coefficients for Ca

coef <- data.frame(Ca=glm_best$coefficients[2,], SV=names(glm_best$coefficients[2,]))

# attaching to the tax table the coef for Calcium. Because there were a lot of rare SVS it does not make sense to visualize the 
# coefficients for those SVs as you cannot really fit a model if you have one observation only. For this reason we will visualize the results 
# for the SVs present in at least 5% of the samples.

ps.sv <- readRDS("../data/Phyloseq_data/ps_sv_PVC.rds")

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


# Gaussian copula models --------------------------------------------------

# Since selection based on AIC is only approximate, we also apply Gaussian Copula models to account for the correlation across taxa. 

mglm_null <- manyglm(sv.MV~ 1, data = samp.data)

mglm_ca <- manyglm(sv.MV~ Ca, data = samp.data, offset = QDoffset)
mglm_ca$AICsum

mglm_status <- manyglm(sv.MV~ Status, data = samp.data, offset = QDoffset)
mglm_status$AICsum

# 53% of models have better AIC compared to Status as covariate
data.frame(ca=mglm_ca$aic, status=mglm_status$aic) %>% 
  mutate(diff=ca-status, res=ifelse(diff<0, "better", "worse"), rows=nrow(.)) %>% 
  filter(res=="better") %>% 
  mutate(rows.filt=nrow(.), perc=rows.filt/rows*100)

# Applying Ecocopula to confirm AIC results
copula.null <- cord(mglm_null)

copula.Ca <- cord(mglm_ca)

copula.Ca$BIC

# factor loading. Read page 396 of Eco-Stats: Data analysis in ecology. Note that these results can change slightly when repeated
ss_ca <- c(sum(copula.null$loadings^2), sum(copula.Ca$loadings^2))

c(ss_ca, 1-ss_ca[2]/ss_ca[1])


## Status
copula.Status <- cord(mglm_status)
copula.Status$BIC

# factor loadings. Read page 396 of Eco-Stats: Data analysis in ecology. Note that these results can change slightly when repeated
ss_st <- c(sum(copula.null$loadings^2), sum(copula.Status$loadings^2))

c(ss_st, 1-ss_st[2]/ss_st[1])

