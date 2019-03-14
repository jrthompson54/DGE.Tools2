library(magrittr)
library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(Xpress2R)


dgeObj <- getRDSobjFromStash("BDL_Rat_LiverSlice_P-20170808-0001_03Dec2017.RDS")

#inspect DGEobj
inventory(dgeObj)
showMeta(dgeObj)

#Get various intensity units
counts <- getItem(dgeObj, "counts")
geneLength <- getItem(dgeObj, "geneData")[["ExonLength"]]
                             
tpm <- convertCounts(counts, unit="tpm", geneLength=geneLength)
log2fpkm <- convertCounts(counts, unit="fpkm", geneLength=geneLength, log=TRUE, normalize="tmm")
log2cpm <- convertCounts(counts, unit="cpm", log=TRUE, normalize="tmm")

#FPR2 data
d_fpr2 <- Xpress2DGEO(21353, level="GRCm38ERCC-ensembl91-genes")


