library(DGEobj)
library(DGE.Tools2)
library(JRTutil)

#omicsoft DGEobj
dge_as <- readRDS("../DGEobj.RDS")
counts <- getItem(dge_as, "counts")
geneData <- getItem(dge_as, "geneData")