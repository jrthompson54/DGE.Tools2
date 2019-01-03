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
dgeObj <- readRDS(file.path(outputPath, "IBD_FLA.RDS"))

counts <- dgeObj$counts_orig
x <- lowIntFilter(counts, zfpkmThreshold = -3, 
                  genelength = dgeObj$geneData_orig$ExonLength,
                  verbose=TRUE)

x <- lowIntFilter(counts, countThreshold = 10, verbose=TRUE)

x <- lowIntFilter(counts, fpkThreshold = 5, 
                  genelength = dgeObj$geneData_orig$ExonLength, verbose=TRUE)


