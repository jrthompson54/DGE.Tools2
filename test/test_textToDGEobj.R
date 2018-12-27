#test CloudToDgeObj function
library(DGE.Tools2)
library(DGEobj)
library(tidyverse)

#generalized path to omicsoft files: need a PID and Omicsoft ProjectName
mountpoint <- "X:"
pid <- "P-20171107-0002"
projectname <- "UCSD_Lung_Fibroblasts_P-20171107-0002_8Feb2018"
path <- file.path(mountpoint, "OmicsoftHome/output", pid, projectname, "ExportedViewsAndTables")
dgeObj <- textToDgeObj(path=path, level="gene", verbose=TRUE)




