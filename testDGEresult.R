RSE <- readRDS("RSE.RDS")
library(magrittr)
library(SummarizedExperiment)
d <- DGEresult()
d %<>% addItem(assay(RSE, "Counts"), "Counts", "assay")
d %<>% addItem(mcols(RSE), "GeneAnnotation", "row")
d %<>% addItem(colData(RSE), "SampAnnotation", "col")

#test trap for overwriting item
d %<>% addItem(colData(RSE), "SampAnnotation", "col")
#force overwrite
d %<>% addItem(colData(RSE), "SampAnnotation", "col", overwrite=T)

#test rmItem
d %<>% rmItem("Counts")

df <- print(d)
