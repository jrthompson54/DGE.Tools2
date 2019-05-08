#### Setup ###
library(stringr)
library(magrittr)
library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(ArrayServer)
library(zFPKM)
library(edgeR)

setwd("~/P01380_BuildDGEobj/IBD_FLA")
inputPath <- "./input"
outputPath <- "./output"
dgeObj <- readRDS(file.path(outputPath, "RNA-Seq_Analysis_of_FL_IBD_Human_Biopsy_P-20170717-0001.RDS"))

dgeObj <- resetDGEobj(dgeObj)


x <- lowIntFilter(dgeObj, zfpkmThreshold = -3,
                  genelength = dgeObj$geneData$ExonLength,
                  verbose=TRUE)

x <- lowIntFilter(dgeObj, countThreshold = 10, verbose=TRUE)

x <- lowIntFilter(dgeObj, fpkThreshold = 5,
                  genelength = dgeObj$geneData$ExonLength, verbose=TRUE)

x <- lowIntFilter(dgeObj, tpmThreshold = 1,
                  genelength = dgeObj$geneData$ExonLength, verbose=TRUE)

x <- lowIntFilter(dgeObj, tpmThreshold = 1, countThreshold=10,
                  genelength = dgeObj$geneData$ExonLength, verbose=TRUE)
