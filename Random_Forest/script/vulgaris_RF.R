library(mlr3)
library(phyloseq)
library(mlbench)
library(tidyverse)
library(metagMisc)
library(Boruta)
library(parallel)
library(rlist)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

getwd()

# Reading file in ---------------------------------------------------------

ps.sv <- readRDS("../../data/Phyloseq_data/ps_sv_PVC.rds")

samp.data <- data.frame(sample_data(ps.sv)) %>% 
  mutate(Sample_ID=str_replace(Sample_ID, "-", ".")) %>% 
  mutate(Dom_event=paste0(.$Dom_event, "_", .$Status))

# make relative abundance and prepare df with subpopulation as target
otu.tab <- data.frame(otu_table(ps.sv)) %>% 
  mutate_if(is.numeric, function(x) {x/sum(x)*100}) %>%
  t() %>%
  data.frame() %>% 
  rownames_to_column(var = "Sample_ID") %>% 
  left_join(., samp.data) %>% 
  select(1:1173, Dom_event, -Sample_ID) %>% 
  mutate(Dom_event=as.factor(Dom_event))

# make relative abundance and prepare df with biological status as target
otu.tab.status <- data.frame(otu_table(ps.sv)) %>% 
  mutate_if(is.numeric, function(x) {x/sum(x)*100}) %>%
  t() %>%
  data.frame() %>% 
  rownames_to_column(var = "Sample_ID") %>% 
  left_join(., samp.data) %>% 
  select(1:1173, Status, -Sample_ID) %>% 
  mutate(Status=as.factor(Status))


# RF Status ---------------------------------------------------------------

boruta.out.status <- Boruta(Status~., otu.tab.status)

# x <- getSelectedAttributes(boruta.out.status, withTentative = F)

# this was saved already
load("../results/Boruta.rda")

# combination of Boruta variables. We are simply sampling combination of features selected by Boruta.

comb <- map(1:5000, ~ map(2:length(x), ~ sample(x,.x)))

out <- list()

# checking combination of SV that could increase the accuracy

for (i in 1:length(comb)) {
  
  out[[i]] <- mclapply(1:length(comb[[1]]), function(x) {
    
    otu.tab.fs <- otu.tab.status %>% 
      select(comb[[i]][[x]], Status)
    
    microbTask <- makeClassifTask(data = otu.tab.fs, target = 'Status')
    rforest <- makeLearner('classif.randomForest')
    getParamSet(rforest)
    dTreeTrained <- train(rforest, microbTask)
    
    kfold_cv <- makeResampleDesc(method = 'RepCV', folds = 5, reps = 10)
    
    
    rforest_cv.status <- mlr::resample(learner = rforest,
                                  task = microbTask,
                                  resampling = kfold_cv,
                                  measures = list(acc, mmce))
    
    
    c(rforest_cv.status$aggr, comb[[i]][[x]])
    
    
  },mc.cores = 8)
  
}

# saveRDS(out, "../results/comb_status.rds")

out <- readRDS("../results/comb_status.rds")

accu <- unlist(map(1:length(out), function(x) {
  map(1:12, function(y) {
    out[[x]][[y]][1]
  })
}))

# maximum accuracy reached
max(accu)
# filtering
out2 <- map(1:length(out), ~ list.filter(out[[.x]], acc.test.mean>=0.93))
# retaining the SVs
indic.spec <- list.clean(out2, function(x) length(x) == 0L, recursive = TRUE)

indic.spec.best <- unlist(indic.spec[[1]])[3:7]

# tuning parameters with all Boruta variables and subselection

vars <- list(x,indic.spec.best)

RF.status <- map(1:2, function(i) {
  
  otu.tab.fs <- otu.tab.status %>% 
    select(all_of(vars[[i]]), Status)
  
  microbTask <- makeClassifTask(data = otu.tab.fs, target = 'Status')
  rforest <- makeLearner('classif.randomForest')
  getParamSet(rforest)
  dTreeTrained <- train(rforest, microbTask)
  
  
  rforestParamSpace <- makeParamSet(
    makeIntegerParam("ntree",lower = 50, upper = 700),
    makeIntegerParam("mtry", lower = 1, upper = 20),
    makeIntegerParam("nodesize", lower = 1, upper = 20)
  )
  
  randSearch <- makeTuneControlRandom(maxit = 100L)
  
  cvForTuning_rf <- makeResampleDesc(method = 'RepCV', 
                                     folds = 5, 
                                     reps = 10)
  
  rforest_tuning.status <- tuneParams(learner = rforest,
                                      task = microbTask,
                                      par.set = rforestParamSpace,
                                      control = randSearch,
                                      resampling = cvForTuning_rf)
  
  
  rforest_tuned.status <- setHyperPars(rforest, par.vals = rforest_tuning.status$x)
  
  kfold_cv <- makeResampleDesc(method = 'RepCV', folds = 5, reps = 10)
  
  
  
  rforest_tuned_cv.status <- mlr::resample(learner =rforest_tuned.status,
                                      task = microbTask,
                                      resampling = kfold_cv,
                                      measures = list(acc, mmce))
  
  
})

# With subselecting at random the variables identified by Boruta, we increased accuracy from 88% to 92%. 
(RF.status[[2]]$aggr[1])-(RF.status[[1]]$aggr[1])


saveRDS(RF.status[[2]], "../results/Status_RF.rds")

Status_RF <- readRDS("../results/Status_RF.rds")

saveRDS(RF.status[[1]], "../results/Status_RF_all.rds")

# plotting the confusion matrix 

Status_RF$pred$data %>% 
  mutate(TR=ifelse(truth==response, TRUE,FALSE), TR=as.factor(TR)) %>% 
  ggplot(aes(truth, TR, color=response)) +
  geom_jitter(alpha=0.7, width = 0.4, height = 15, size=1)+
  ylab("")+
  theme_Publication()+
  scale_color_manual(values = c("Red3", "darkgreen"))+
  theme(legend.position = "none", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  coord_flip()


# Now we repeat the same for the domestication event within domestication status. Basically, AE, AD, MD, MW

# RF Dom_event ------------------------------------------------------------

boruta.out <- Boruta(Dom_event~., otu.tab)

x <- getSelectedAttributes(boruta.out, withTentative = F)

# save(x, file = "../results/Boruta_dom_event.rda")

load("../results/Boruta_dom_event.rda")

# combination of Boruta variables

comb <- map(1:5000, ~ map(2:length(x), ~ sample(x,.x)))

out <- list()

for (i in 1:length(comb)) {
  
  out[[i]] <- mclapply(1:length(comb[[1]]), function(x) {
    
    otu.tab.fs <- otu.tab %>% 
      select(comb[[i]][[x]], Dom_event)
    
    microbTask <- makeClassifTask(data = otu.tab.fs, target = 'Dom_event')
    rforest <- makeLearner('classif.randomForest')
    dTreeTrained <- train(rforest, microbTask)
    
    kfold_cv <- makeResampleDesc(method = 'RepCV', folds = 5, reps = 10)
    
    
    rforest_cv.status <- mlr::resample(learner = rforest,
                                  task = microbTask,
                                  resampling = kfold_cv,
                                  measures = list(acc, mmce))
    
    
    c(rforest_cv.status$aggr, comb[[i]][[x]])
    
    
  },mc.cores = 8)
}

saveRDS(out, "../results/comb_dom_event.rds")

out <- readRDS("../results/comb_dom_event.rds")

accu <- unlist(map(1:length(out), function(x) {
  map(1:13, function(y) {
    out[[x]][[y]][1]
  })
}))

# max accuracy reached is 0.75

max(accu)

out2 <- map(1:length(out), ~ list.filter(out[[.x]], acc.test.mean>=0.75))

indic.spec <- list.clean(out2, function(x) length(x) == 0L, recursive = TRUE)

indic.spec.best <- unlist(indic.spec[[1]])[3:13]

# tuning parameters with all Boruta variables and subselection

vars <- list(x,indic.spec.best)

RF.dom_event <- map(1:2, function(i) {
  
  otu.tab.fs <- otu.tab %>% 
    select(all_of(vars[[i]]), Dom_event)
  
  microbTask <- makeClassifTask(data = otu.tab.fs, target = 'Dom_event')
  rforest <- makeLearner('classif.randomForest')
  getParamSet(rforest)
  dTreeTrained <- train(rforest, microbTask)
  
  
  rforestParamSpace <- makeParamSet(
    makeIntegerParam("ntree",lower = 50, upper = 700),
    makeIntegerParam("mtry", lower = 1, upper = 20),
    makeIntegerParam("nodesize", lower = 1, upper = 20)
  )
  
  randSearch <- makeTuneControlRandom(maxit = 100L)
  
  cvForTuning_rf <- makeResampleDesc(method = 'RepCV', 
                                     folds = 5, 
                                     reps = 10)
  
  rforest_tuning.status <- tuneParams(learner = rforest,
                                      task = microbTask,
                                      par.set = rforestParamSpace,
                                      control = randSearch,
                                      resampling = cvForTuning_rf)
  
  
  rforest_tuned.status <- setHyperPars(rforest, par.vals = rforest_tuning.status$x)
  
  kfold_cv <- makeResampleDesc(method = 'RepCV', folds = 5, reps = 10)
  
  
  rforest_tuned_cv.status <- mlr::resample(learner =rforest_tuned.status,
                                      task = microbTask,
                                      resampling = kfold_cv,
                                      measures = list(acc, mmce))
  
  
})


(RF.dom_event[[2]]$aggr[1])-(RF.dom_event[[1]]$aggr[1])

saveRDS(RF.dom_event[[2]], "../results/Dom_event_RF.rds")

Dom_event_RF <- readRDS("../results/Dom_event_RF.rds")

saveRDS(RF.dom_event[[1]], "../results/Dom_event_RF_all.rds")


# plotting the confusion matrix 

Dom_event_RF$pred$data %>%
  mutate(truth=case_when(truth=="Mesoamerica_Wild" ~ "MW",
                         truth=="Mesoamerica_Domesticated" ~ "MD", 
                         truth=="Ande_Domesticated" ~ "AD",
                         truth=="Ande_Wild" ~ "AW"), response=case_when(response=="Mesoamerica_Wild" ~ "MW",
                                                                      response=="Mesoamerica_Domesticated" ~ "MD", 
                                                                      response=="Ande_Domesticated" ~ "AD",
                                                                      response=="Ande_Wild" ~ "AW")) %>% 
  mutate(truth=fct_relevel(truth, "AW", "MW", "AD", "MD")) %>% 
  mutate(TR=ifelse(truth==response, TRUE,FALSE), TR=as.factor(TR)) %>% 
  ggplot(aes(truth, TR, color=response)) +
  geom_jitter(alpha=0.7, width = 0.4, height = 15, size=1)+
  ylab("")+
  theme_Publication()+
  theme(legend.position = "none", 
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  coord_flip()
