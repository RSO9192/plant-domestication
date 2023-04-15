# This script was used to calculate the difference in predicted functional profiles between wild and domesticated accessions for the 2nd experiment (P. vulgaris grown in the field). 
# The same logic applies for the other experiment. The in-house genomic database can be requested

library(Tax4Fun2)
library(tidyverse)
library(metagMisc)
library(phyloseq)
library(vegan)

# prepare metadata
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

# otu table with no contaminant, only bacterial reads
ps <- readRDS("../data/Phyloseq_data/ps_sv_PVC.rds")
ps <- phyloseq_filter_prevalence(ps, prev.trh = 0.01, abund.trh = 20, threshold_condition = "OR", abund.type = "total")

otu.table <- as.data.frame(otu_table(ps)) %>% 
  add_rownames(var="seq") %>% 
  mutate(SV=paste0("SV_", 1:n()))

metadata <- data.frame(sample_data(ps)) %>% 
  rename(Accession=Sample_ID)

rownames(metadata) <- NULL

# adonis ------------------------------------------------------------------
adon <- read.delim("../data/functions_data/functional_prediction_PVC.txt", sep = "\t") %>% 
  select(-description) %>% 
  t()

colnames(adon) <- adon[1,]

adon <- as.data.frame(adon[-1,]) %>% 
  mutate_all(as.numeric)

adonis2(adon ~ Status, data=metadata, method = "robust.aitchison") # used because data is compositional
