#test runVoom
rm(list=ls())
# dgelist <- readRDS('~/R/recount/SRP010041/dgelist.RDS')
# load(file.path('~/R/recount/SRP010041', 'rse_gene.Rdata'))
library(dplyr)
library(edgeR)
library(limma)
library(magrittr)
library(DGEobj)
library(DGE.Tools2)

setwd("~/R/DGE.Tools_Example")

dgeObj <- OmicsoftToDgeObj()

dgeObj <- runEdgeRNorm(dgeObj)

#define a formula and construct a design matrix
design <- getItem(dgeObj, "design")
design$ReplicateGroup %<>% as.factor
design$ReplicateGroup %<>% relevel("Normal_control")
formula <- '~ ReplicateGroup'
designMatrix <- model.matrix (as.formula(formula), design)

#try all 6 senarios
d1 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = FALSE)

#QW and Var.design
vd <- model.matrix(as.formula("~ Treatment"), design)
d2 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = FALSE,
              var.design=vd)

#QW and Var.design and dupCor
block <- c(1,2,3,1,2,3,4,5,6,4,5,6,7,8,9,7,8,9)
d3 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = FALSE,
              var.design=vd,
              dupcorBlock=block)
#add QW
d4 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = TRUE)

#QW and Var.design
vd <- model.matrix(as.formula("~ Treatment"), design)
d5 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = TRUE,
              var.design=vd)

#QW and Var.design and dupCor
block <- c(1,2,3,1,2,3,4,5,6,4,5,6,7,8,9,7,8,9)
d6 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = TRUE,
              var.design=vd,
              dupcorBlock=block)


