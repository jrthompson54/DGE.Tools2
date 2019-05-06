library(magrittr)
library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)

dgeObj <- readRDS("Z:/DGEobj_library(stash copy)/ColumbiaHearts_P-20151002-0004_R90_27Nov2017.RDS")

PreVsRef <- dgeObj$PreVsRef
PostVsRef <- dgeObj$PostVsRef


ttList <- getType(dgeObj, "topTable")[1:2]
names(ttList)

cPlotDat <- comparePrep(ttList)

