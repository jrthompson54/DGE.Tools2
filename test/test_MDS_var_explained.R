
library(magrittr)
library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(limma)

mdsresult <- readRDS("Z:/SomaDat/output/mdsPlot_x079.RDS")[[2]]

varResults <- MDS_var_explained(mdsresult, baseFontSize = 14, barWidth=0.6)
varResults[[1]] #the Variance per dimension plot
varResults[[2]] #the cumulative variance plot
var_explained <- varResults[[3]]  #data used for plotting (unfiltered)


setBreaks <- function(limits){
  #return integer breaks
  low <- floor(limits[1])
  high <- ceiling(limits[2])
  seq(from=low, to=high, by=1)
}

varResults[[1]] + scale_x_continuous(breaks = setBreaks)

######################################################
