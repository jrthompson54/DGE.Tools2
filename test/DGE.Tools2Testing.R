#DGE.Tools2 testing
#
library(gdata)
library(openxlsx)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
library(reshape2)
library(ggplot2)
library(assertthat)
library(magrittr)
library(DGEobj)
library(DGE.Tools2)

setwd("~/R/DGE.Tools_Example")

#test OmicsoftToDgeObj
d <- OmicsoftToDgeObj(customAttr = list(
    Genome="Human.B38",
    GeneModel="Ensembl.R82"))
class(d)
names(d)

#test Build_RSE
rse <- Build_RSE()

#test convertCounts
mycounts <- getItem(d, "counts")





    
