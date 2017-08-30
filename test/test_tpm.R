library(DGEobj)
library(DGE.Tools2)
library(Xpress2R)
library(magrittr)
library(tidyverse)
#Get CCL4 data from Omicsoft
do <- readRDS("~\\R\\contrastDB\\DGEobj\\CCL4_DGEobj.RDS")
do <- annotateDGEobj(do, "~\\R\\contrastDB\\DGEobj\\CCL4_P-20160407-0002_20Jul2016_FastqRegistration.txt")
attr(do, "PID") 
names(attributes(do))
dim(do)

#Get XpressData
if (file.exists("X19955.RDS")) {
    dx <- readRDS("X19955.RDS")
} else {
    dx <- Xpress2DGEO(19955)
    saveRDS (dx, "X19955.RDS")
}

#print length==1 attributes in a table
x <- attributes(dx)
idx <- lapply(x, length) == 1
y <- stack(x[idx])
colnames(y) <- c("Value", "Attribute")
y <- select(y, Attribute, Value)


counts.dx <- getItem(dx, "counts_orig")
counts.do <- getItem(do, "counts_orig")
ga.dx <- getItem(dx, "geneData_orig")
ga.do <- getItem(do, "geneData_orig")

#should have same dim
tpm.do <- tpm(do)
tpm2.do <- convertCounts(counts.do, unit="tpm", geneLength=ga.do$ExonLength)
dim(tpm.do)
dim(counts.do)
inspect(counts.do)
inspect(tpm.do)
#looks consistent.  Seems like a high degree of zeros (~70%).  Maybe this is a mRNA dataset???




#now try with Xpress
tpm.dx <- tpm(dx)
dim(tpm.dx)
dim(counts.dx)
inspect(counts.dx)
inspect(tpm.dx)
