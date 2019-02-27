# voomWorkflow
library(magrittr)
library(tidyverse)
library(DGEobj)
library(JRTutil)


 #as it appears in Omicsoft and the Regfile
#Paths to Omicsoft files follow a convention requiring the PID and the Omicsoft Project name.
mountpoint <- "y:"
pid <- P-20161026-0001
projectName <- "MOGAT2_Inhib_MCD-HFD_P-20161026-0001_29Nov2016"
dgeObj <- buildO


dgeObj <- textToDgeObj(path=path, level="gene", verbose=TRUE)
projectName <- "RNA-Seq_Analysis_of_FL_IBD_Human_Biopsy_P-20170717-0001" #as it appears in Omicsoft and the Regfile

rdsname <- "RNA-Seq_Analysis_of_FL_IBD_Human_Biopsy_P-20170717-0001"
regfile <- file.path(inputPath, str_c(projectName, ".txt"))  #omicsoft registration file; .txt or .gz file
designMatrixName <- "Category2"   #appropriate name for this model
formula <- '~ 0 + Category2'


dgeObj <- getRDSobjFromStash("MOGAT2_Inhib_MCD-HFD_P-20161026-0001_29Nov2016.RDS")
dgeObj <- getRDSobjFromStash(stringr::str_c(projectName, ".RDS"))
#Use Design$TRDSample.x for this dataset
dupcorBlock <- NULL    #define duplicates for dupliceCorrelation method; set to NULL to disable

#strip down DGEobj
#d <- dgeObj
