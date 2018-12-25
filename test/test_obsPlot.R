# Test obsPlot2.R


library(DGEobj)
library(DGE.Tools2)
library(tidyverse)
library(magrittr)
library(JRTutil)

stashRoot <- getStashLocation()
stashPath <- file.path(stashRoot, "data/nonclin/DGEobj_library")
dgeObj <- readRDS(file.path(stashPath, "BDL_Rat_LiverSlice_P-20170808-0001_03Dec2017.RDS"))

tcDat <-tidyContrasts(dgeObj, rownameColumn="EnsgID")
unique(tcDat$contrast)
#now plot logration +/- CI


#Intensity Plots
log2cpm <- convertCounts(dgeObj$counts, unit="cpm", log=TRUE, normalize="tmm")

#get EnsgID to GeneSymbol map
ga <- dgeObj$geneData %>%
  rownames_to_column(var="EnsgID") %>%
  select(EnsgID, GeneSymbol=GeneName)

#filter for genes of interest
##let's plot lpar gene family
idx <- str_detect(ga$GeneSymbol, "^Lpar")
GOI <- ga$GeneSymbol[idx]
log2cpm <- log2cpm[idx,]

tiDat <- tidyIntensity(log2cpm, rowIDcolumn="EnsgID", dgeObj$design$ReplicateGroup)
#add gene symbols
tiDat %<>% left_join(ga)


myplot <- obsPlot2(tiDat, plotByCol="GeneSymbol",
                   groupCol="group",
                   valueCol="log2cpm",
                   scale="fixed")
myplot
