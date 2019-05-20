# voomWorkflow
library(magrittr)
library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(Xpress2R)


#Test1 Omicsoft project
#
#projectName as it appears in Omicsoft and the Regfile
# projectName <- "omicsoft/TB3_TGFB_Stim_HuIPF_Fibroblasts_Ensembl_30Jun2018"
# dgeObj <- buildOmicsoftDGEobj(projectName, mountPoint="y:", PID="TB3", omicsoftHomePath = "omicsoft-legacy/omicsoft_ngs")
mountPoint <- "y:"
projectName <- "BDL_Rat_LiverSlice_P-20170808-0001_03Dec2017"
dgeObj <- JRTutil::buildOmicsoftDGEobj(projectName = projectName,
                                       level="gene",
                                       mountPoint=mountPoint)
dgeObj$design$ReplicateGroup <- dgeObj$design$replicate_group

dim(dgeObj)
inventory(dgeObj)

designMatrixName <- "ReplicateGroup"   #appropriate name for this model
formula <- "~ 0 + ReplicateGroup"

dupcorBlock <- NULL    #define duplicates for dupliceCorrelation method; set to NULL to disable

d <- voomWorkflow(dgeObj,
                  formula=formula,
                  projectName = projectName,
                  designMatrixName = designMatrixName,
                  outputPath = "z:/testOutput",
                  countThreshold = 10,
                  tpmThreshold = 1)
dim(d)
inventory(d)



#test2 Xpress
d2 <- getRDSobjFromStash("RNA-Seq_Analysis_of_FL_IBD_Human_Biopsy_P-20170717-0001.RDS")
#Use Design$TRDSample.x for this dataset
dupcorBlock <- NULL    #define duplicates for dupliceCorrelation method; set to NULL to disable
