#isoformFrac test

library(DGEobj)
library(DGE.Tools2)
setwd("C:/Users/thompj27/Documents/Fibrosis/FGF21/Stelic FGF21/3rdStudy/Transcript")

d <- OmicsoftToDgeObj(level="isoform")

newd <- isoformFrac(d)
