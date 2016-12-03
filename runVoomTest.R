#test runVoom

dgelist <- readRDS('~/R/recount/SRP010041/dgelist.RDS')
load(file.path('~/R/recount/SRP010041', 'rse_gene.Rdata'))
library(dplyr)
library(edgeR)
library(limia)
source('runVoom.R')

design <- colData(rse_gene) %>% as.data.frame
design$State <- as.factor(c(rep("Norm", 3), rep("IPF", 3)))
design$State %<>% relevel("Norm")
formula <- '~ State'
designMatrix <- model.matrix (as.formula(formula), design)


dgeResult <- runVoom(dgelist, designMatrix, formula, qualityWeights = FALSE)

result <- runVoom(dgelist, designMatrix, formula,
                    dupcorBlock=NULL,
                    qualityWeights = TRUE,
                    qwBlock=NULL,
                    resultList=NULL)
