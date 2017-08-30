#Somethings wierd with the Xpress counts.
#Swap in a Omicsoft Human Ensembl dataset for the test

#tpm function testing

library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(magrittr)
library(edgeR)

#Xpress DGEO
dx <- readRDS("~\\IO\\H&N\\RNA-Seq_P-20170209-0006\\input\\X20213_Gene_DGEO.RDS")
#Omicosft DGEO
do <- readRDS("~\\Fibrosis\\P00785_Metalloprotease_Cell_Type_Analysis\\input\\LungFib_DGEobj.RDS")
attr(do, "source") <- "Omicsoft"
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

