library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(tidyverse)
library(magrittr)

dgeObj <- readRDS("./output/IPF_Lung_UPenn_P-20161128-0001_11Sep2017_DGEobj.RDS")
dim(dgeObj)

newDgeObj <- resetDGEobj(dgeObj)
dim(newDgeObj)
inventory(newDgeObj)
