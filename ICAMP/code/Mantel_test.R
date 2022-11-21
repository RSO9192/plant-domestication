library(vegan)
library(phyloseq)
library(metagMisc)
library(tidyverse)
library(readxl)
source("tools.r")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

getwd()

ps.sv <- readRDS("../../3.filtering/results/tables/ps.rds")

# Removing SVs that are present at 1% and have total count lower than 20
ps.sv <- phyloseq_filter_prevalence(ps.sv, prev.trh = 0.01, abund.trh = 20, threshold_condition = "OR", abund.type = "total")

# inclusind collection time
collection_time <- read_excel("../../../../../../../Accessions to choose/datasets/Stations_collection_year.xlsx") %>% 
  select(ACCESION, 'COLLECTED YEAR') %>% 
  rename(Sample_ID=ACCESION, C.year='COLLECTED YEAR') %>% 
  mutate(Sample_ID=str_squish(Sample_ID)) %>% 
  mutate(Sample_ID=str_remove(Sample_ID, " "), C.year=as.numeric(C.year)) %>% 
  mutate(C.year.d=C.year-2019) %>% 
  mutate(C.year.d=C.year.d*-1)

samp.data <- data.frame(sample_data(ps.sv))

samp.data <- left_join(samp.data, collection_time) %>% 
  select(Sample_ID,C.year,Ca,Altitude, Flowering_time, Latitude, Longitude,12:24)

rownames(samp.data) <- samp.data$Sample_ID

samp.data <- samp.data[,-1, drop=FALSE] 

#mutate(C.year=log(C.year))

# environmental vector to test- euc distances

icres <- readRDS("../outputs/ICAMP.detail.rda")
icres <- icres[[1]]

time=(samp.data[match(icres$sample1,rownames(samp.data)),]+
        samp.data[match(icres$sample2,rownames(samp.data)),])/2

#colnames(menv)=paste0("m",colnames(env.use))


denv=abs(samp.data[match(icres$sample1,rownames(samp.data)),]-
           samp.data[match(icres$sample2,rownames(samp.data)),])

mdenv=cbind(icres[,1:2],denv)

# to account for sampling and location

grp.rand=data.frame(sample_data(ps.sv)) %>% 
left_join(., collection_time) %>% 
  select(C.year, Regeneration_site, Sample_ID) %>% 
  mutate(var=paste0(C.year, "_", Regeneration_site)) %>% 
  select(var, Sample_ID) %>% 
column_to_rownames("Sample_ID")

unique(grp.rand$var)

# read in pair wise comparison and perform mantel test

out <- map(3:7, function(i) {
  df <- icres %>% 
    select(1,2,i) # HoS
  
  map(colnames(mdenv)[4:18], function(x) {
    
    mdenv2 <- mdenv %>% 
      select(1:2, x)
    
    mci=mcMantel(y.3col = df,y.ids = NULL,
                 x.3col = mdenv2, grp.rand = NULL,
                 grp.const = grp.rand,try.time = 4, na.rm = TRUE,
                 method = "spearman",permutations = 999)  
    
    c(mci$statistic, mci$signif)
    
})})

names(out) <- colnames(icres)[3:7]

for (i in 1:length(out)) {
  
names(out[[i]]) <- colnames(mdenv[4:18])
  
}

saveRDS(out, "../results/out_Mantel.rds")

out <- readRDS("../results/out_Mantel.rds")

out.df <- data.frame(Vb=colnames(mdenv)[4:18]) %>% 
  add_column(HeS=NA, HoS=NA, DL=NA, DH=NA, DR=NA)

out.df[3,2] <- -0.07
out.df[5,4] <- 0.18
out.df[8,4] <- -0.14
out.df[6,5] <- -0.12
out.df[5,5] <- 0.13


