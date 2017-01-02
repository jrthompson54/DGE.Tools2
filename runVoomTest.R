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

#QW and Var.design with eBayes
vd <- model.matrix(as.formula("~ Treatment"), design)
d7 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = TRUE,
              var.design=vd, 
              runEBayes=TRUE)

m <- ggplotMDS(dgeObj, colorBy = design$ReplicateGroup, labels=design$ReplicateGroup)

m <- ggplotMDS(dgeObj, colorBy = design$ReplicateGroup, symSize =5)
m[[1]]

#use color and shape
m <- ggplotMDS(dgeObj, colorBy = design$ReplicateGroup, 
               shapeBy = design$Disease.Status, symSize =5)
m[[1]]

#use color, shape and size
m <- ggplotMDS(dgeObj, colorBy = design$ReplicateGroup, 
               shapeBy = design$Disease.Status, 
               sizeBy = design$Treatment)
m[[1]]

#use color and shape with labels and labelSize
m <- ggplotMDS(dgeObj, colorBy = design$Treatment, 
               shapeBy = design$Disease.Status, symSize =5,
               labels=design$VendorBarcode,
               labelSize=3)
m[[1]]

#use color and shape without labels
m <- ggplotMDS(dgeObj, colorBy = design$Treatment, 
               shapeBy = design$Disease.Status, symSize =5,
               labels=NULL)

