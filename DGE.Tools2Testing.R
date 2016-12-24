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
library(DGEobj)
library(magrittr)
library(DGE.Tools2)

#test OmicsoftToDgeObj
#
setwd("~/R/DGE.Tools_Example")

d <- OmicsoftToDgeObj(customAttr = list(
    Genome="Human.B38",
    GeneModel="Ensembl.R82"))
    
