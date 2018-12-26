# Test tidyContrasts, logRatioPlot, tidyIntensity and obsPlot2


library(DGEobj)
library(DGE.Tools2)
library(tidyverse)
library(magrittr)
library(JRTutil)

stashRoot <- getStashPath()
stashRoot <- "//stash.pri.bms.com/stash"
stashPath <- file.path(stashRoot, "data/nonclin/DGEobj_library")
dgeObj <- readRDS(file.path(stashPath, "BDL_Rat_LiverSlice_P-20170808-0001_03Dec2017.RDS"))

#get EnsgID to GeneSymbol map
ga <- dgeObj$geneData %>%
  rownames_to_column(var="EnsgID") %>%
  select(EnsgID, GeneSymbol=GeneName)

tcDat <-tidyContrasts(dgeObj, rownameColumn="EnsgID", includeColumns = c("logFC", "CI.R", "CI.L"))
#add gene symbols to tcDat
tcDat %<>% left_join(ga)

##let's plot lpar gene family
idx <- str_detect(ga$GeneSymbol, "^Lpar")
GOI <- ga$GeneSymbol[idx]
tcDat %<>% filter(GeneSymbol %in% GOI)
unique(tcDat$Contrast)

#now plot logratio +/- CI
MyPlot <- logRatioPlot(tcDat, plotType="point",
                       facetColname = "GeneSymbol",
                       xColname = "Contrast",
                       facetCol=4,
                       scales="fixed",
                       facet=TRUE,
                       title = "Test",
                       pointSize=4,
                       lineLayer=TRUE,
                       lineSize=0.1,
                       xAngle=60,
                       ylab = "Log2Ratio")

#now plot logratio +/- CI
MyPlot <- logRatioPlot(tcDat,
                       facetColname = "GeneSymbol",
                       xColname = "Contrast")

#Intensity Plots
log2cpm <- convertCounts(dgeObj$counts, unit="cpm", log=TRUE, normalize="tmm")



#filter for genes of interest
##let's plot lpar gene family
idx <- str_detect(ga$GeneSymbol, "^Lpar")
GOI <- ga$GeneSymbol[idx]
log2cpm <- log2cpm[idx,]

tiDat <- tidyIntensity(log2cpm,
                       rowidColname ="EnsgID",
                       keyColname="Contrast",
                       valueColname = "Log2CPM",
                       group=dgeObj$design$ReplicateGroup)
#add gene symbols
tiDat %<>% left_join(ga)


myplot <- obsPlot2(tiDat, plotByCol="GeneSymbol",
                   groupCol="group",
                   valueCol="Log2CPM",
                   scale="fixed")
myplot
