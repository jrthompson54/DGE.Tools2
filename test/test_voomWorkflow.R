# voomWorkflow
library(magrittr)
library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)



#projectName as it appears in Omicsoft and the Regfile
# projectName <- "MOGAT2_Inhib_MCD-HFD_P-20161026-0001_29Nov2016"
# dgeObj <- buildOmicsoftDGEobj(projectName, mountPoint="y:")


projectName <- "RNA-Seq_Analysis_of_FL_IBD_Human_Biopsy_P-20170717-0001"
# regfile <- file.path(inputPath, str_c(projectName, ".txt"))  #omicsoft registration file; .txt or .gz file
designMatrixName <- "Category2"   #appropriate name for this model
formula <- '~ 0 + Category2'

#Use Design$TRDSample.x for this dataset
dupcorBlock <- NULL    #define duplicates for dupliceCorrelation method; set to NULL to disable

dgeObj <- getRDSobjFromStash(stringr::str_c(projectName, ".RDS"))
dgeObj <- voomWorkflow(dgeObj,
                       formula=formula,
                       projectName = projectName,
                       designMatrixName = designMatrixName,
                       outputPath = "z:/testOutput")
