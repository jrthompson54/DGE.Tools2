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

setwd("~/R/lib/pkgsrc/DGE.Tools2/inst/extdata")
setwd("~/R/lib/pkgsrc/DGE.Tools2")
dgeObj <- OmicsoftToDgeObj(path="./inst/extdata")
# dgeObj <- OmicsoftToDgeObj()

setwd("~/R/lib/pkgsrc/DGE.Tools2")
dgeObj <- runEdgeRNorm(dgeObj)

#low expression filter
counts <- getItem(dgeObj, "counts")
genelength <-getItem(dgeObj, "geneData")$ExonLength
fpk <- convertCounts(counts, 
                     unit="FPK", 
                     geneLength=genelength)
#keep FPK >=5 in 50% of samples
idx <- fpk >= 5.0
frac <- rowSums(fpk) / ncol(fpk)
idx <- frac >= 0.75
d1 <- dgeObj  #save for testing
dgeObj <- subset(dgeObj, idx, 1:ncol(dgeObj))
d2 <- dgeObj

#define a formula and construct a design matrix
design <- getItem(dgeObj, "design")
design$ReplicateGroup %<>% as.factor
design$ReplicateGroup %<>% relevel("Normal_control")
formula <- '~ 0 + ReplicateGroup'
#build the designMatrix and add some attributes
designMatrix <- model.matrix (as.formula(formula), design)
designMatrixName <- "Treatment"
#add some attributes
designMatrix <- setAttributes(designMatrix, list(formula=formula, 
                                 parent="design"))
#save the designMatrix
dgeObj <- addItem(dgeObj, designMatrix, designMatrixName, "designMatrix")


#dispersion plot
dispPlot <- plotDisp(getItem(dgeObj, "DGEList"), designMatrix)
dispPlot + baseTheme(18)


#QW and Var.design and dupCor
block <- c(1,2,3,1,2,3,4,5,6,4,5,6,7,8,9,7,8,9)
vd <- model.matrix(as.formula("~ Treatment"), design)
d1 <- runVoom(dgeObj, designMatrixName, 
              qualityWeights = TRUE,
              var.design=vd,
              dupcorBlock=block)

#QW and Var.design
# vd <- model.matrix(as.formula("~ Treatment"), design)
# d2 <- runVoom(dgeObj, designMatrixName, qualityWeights = FALSE,
#               var.design=vd)

#use color and shape with labels and labelSize
m <- ggplotMDS(d1, colorBy = design$Treatment, 
               shapeBy = design$Disease.Status, symSize =5,
               labels=design$VendorBarcode,
               labelSize=3)
m[[1]]

# source('C:/Users/thompj27/Documents/R/lib/pkgsrc/DGE.Tools2/R/runContrasts.R')
#runContrast testing  
contrastList  <- list(TGF_Norm = "ReplicateGroupNormal_TGFb - ReplicateGroupNormal_control",
                      TGF_Stable = "ReplicateGroupStable_TGFb - ReplicateGroupStable_control",
                      TGF_Rapid = "ReplicateGroupRapid_TGFb - ReplicateGroupRapid_control"
)
# library(assertthat)
DgeObj_contrast <- runContrasts(d1, "Treatment_fit", contrastList, runTopTreat=T)
saveRDS(DgeObj_contrast, "../DGEobj.RDS")

#sva test
dgeObj_sva <- runSVA(d1, "Treatment")


#Qvalue test
#extract a topTable List
contrastList <- getType(DgeObj_contrast, type="topTAble", parent="Treatment_fit_cf")
contrastList2 <- runQvalue(contrastList)
#now put the contrast back in the DGEobj
DgeObj_q <- addItems(DgeObj_contrast, contrastList2, overwrite=TRUE)

