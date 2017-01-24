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

setwd("~/R/lib/pkgsrc/DGE.Tools2")

dgeObj <- OmicsoftToDgeObj(path="./inst/extdata")

dgeObj <- runEdgeRNorm(dgeObj)

#define a formula and construct a design matrix
design <- getItem(dgeObj, "design")
design$ReplicateGroup %<>% as.factor
design$ReplicateGroup %<>% relevel("Normal_control")
formula <- '~ 0 + ReplicateGroup'
# designMatrix <- model.matrix (as.formula(formula), design)

#try 2 senarios
#Simple Voom (SV)
d1 <- runVoom(dgeObj, formula, 
              formulaName = "SV",
              qualityWeights = FALSE)

#Var.design (SV_VD)
vd <- model.matrix(as.formula("~ Treatment"), design)
d2 <- runVoom(dgeObj, formula, 
              formulaName = "SV_VD",
              qualityWeights = FALSE,
              var.design=vd)

#Var.design and dupCor
block <- c(1,2,3,1,2,3,4,5,6,4,5,6,7,8,9,7,8,9)
d3 <- runVoom(dgeObj, formula, 
              formulaName = "SV_VD_DC",
              qualityWeights = FALSE,
              var.design=vd,
              dupcorBlock=block)
#add QW
d4 <- runVoom(dgeObj, formula, 
              formulaName = "VQW",
              qualityWeights = TRUE)

#QW and Var.design
vd <- model.matrix(as.formula("~ Treatment"), design)
d5 <- runVoom(dgeObj, formula, 
              formulaName = "VQW_VD",
              qualityWeights = TRUE,
              var.design=vd)

#QW and Var.design and dupCor
block <- c(1,2,3,1,2,3,4,5,6,4,5,6,7,8,9,7,8,9)
d6 <- runVoom(dgeObj, formula, 
              formulaName = "VQW_VD_DC",
              qualityWeights = TRUE,
              var.design=vd,
              dupcorBlock=block)

#QW and Var.design with eBayes
block <- c(1,2,3,1,2,3,4,5,6,4,5,6,7,8,9,7,8,9)
vd <- model.matrix(as.formula("~ Treatment"), design)
d7 <- runVoom(dgeObj, formula, 
              formulaName = "VQW_VD_DC",
              qualityWeights = TRUE,
              var.design=vd,
              dupcorBlock=block,
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
m[[1]]

#use color and shape default labels (should be ReplicateGroup
m <- ggplotMDS(dgeObj, colorBy = design$Treatment, 
               shapeBy = design$Disease.Status, symSize =5)
m[[1]]


#runContrast testing  
contrastList  <- list(TGF_Norm = "ReplicateGroupNormal_TGFb - ReplicateGroupNormal_control",
                      TGF_Stable = "ReplicateGroupStable_TGFb - ReplicateGroupStable_control",
                      TGF_Rapid = "ReplicateGroupRapid_TGFb - ReplicateGroupRapid_control"
)
library(assertthat)
DgeObj_contrast <- runContrasts(d7, contrastList)
