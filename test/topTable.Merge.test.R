library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(tidyverse)
d <- readRDS("./BDL_DGEobj.RDS")

ttl <- getType(d, "topTable")

x <- extractCol(ttl, "P.Value")
y <- topTable.merge(ttl)
