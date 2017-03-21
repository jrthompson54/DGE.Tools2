#test convertCounts with vector

library(DGEobj)
library(DGE.Tools2)

#Load DGE.Tools2 sample data
rawDataPath <- paste(.libPaths()[[1]], "/DGE.Tools2/extdata", sep="")
d <- OmicsoftToDgeObj(path=rawDataPath)

counts <- getItem(d, "counts")
#get a numeric vector
c1 <- counts[,1]
class(c1)
ga <- getItem(d, "geneData")

fpkm <- convertCounts(c1, unit = "FPKM", geneLength = ga$ExonLength)
