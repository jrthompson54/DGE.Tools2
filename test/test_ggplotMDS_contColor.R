#  test ggplotMDS with continuous color scale
library(magrittr)
library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)

dgeObj <- getRDSobjFromStash("BDL_Rat_LiverSlice_P-20170808-0001_03Dec2017.RDS")
dgeObj$design$number <- 1:ncol(dgeObj)
log2CPM <- convertCounts(dgeObj$counts, unit="cpm", log=TRUE, normalize = "tmm")
mdsresult <- ggplotMDS(log2CPM, colorBy = dgeObj$design$number, shapeBy=dgeObj$design$compound, labels=NULL)
mdsresult$plot

