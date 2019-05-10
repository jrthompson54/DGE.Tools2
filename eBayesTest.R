### Autoprop testing

#### Setup ###
library(stringr)
library(magrittr)
library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(ArrayServer)
library(zFPKM)
library(edgeR)

if (!file.exists("~/eBayesTest")) dir.create("~/eBayesTest")
setwd("~/eBayesTest")
inputPath <- "./input"
outputPath <- "./output"
dgeObj <- getRDSobjFromStash("MOGAT2_Inhib_MCD-HFD_P-20161026-0001_29Nov2016.RDS")

design <- dgeObj$design
cfit <- dgeObj$Category2_fit_cf
ttlist <- getType(dgeObj, "topTable")
cat("Standard pipeline prop=0.01\n")
knitr::kable(summarizeSigCounts(ttlist))

#run contrasts with different prop setting (try 500, 5000 of 20000 i.e. 0.025 and 0.25)
############################
## Rerun with proportion = 0.01
############################

d2 <- rmItems(dgeObj, items=c(11:16))


d_prop01 <- runVoom(d2,
                    designMatrixName="ReplicateGroup",
                    runEBayes = FALSE,
                    mvPlot=FALSE)

contrastList  <- list(CDAHFD_vs_Control = "ReplicateGroupMCD.HFD_Veh - ReplicateGroupControl_Veh",
                      BMS963272_3_vs_CDAHFD = "ReplicateGroupMCD.HFD_BMS963272_Hi - ReplicateGroupMCD.HFD_Veh"
)
d_prop01 <- runContrasts(d_prop01,
                         designMatrixName="ReplicateGroup",
                         contrastList=contrastList,
                         runTopTreat=TRUE,
                         FoldChangeThreshold = 2,
                         Qvalue=TRUE,
                         IHW=TRUE,
                         runEBayes=TRUE,
                         proportion=0.01,
                         verbose=TRUE
                         )

cat("Proportion = 0.01 Results:\n")
cat("No FC Threshold\n")
knitr::kable(summarizeSigCounts(getType(d_prop01, "topTable")))
cat("FC Threshold = 2\n")
knitr::kable(summarizeSigCounts(getType(d_prop01, "topTable"), fcThreshold = 2))
cat("treat FC Threshold = 2\n")
knitr::kable(summarizeSigCounts(getType(d_prop01, "topTreat")))

############################
## Rerun with proportion = 0.05
############################

d2 <- rmItems(dgeObj, items=c(11:16))


d_prop05 <- runVoom(d2,
              designMatrixName="ReplicateGroup",
              runEBayes = FALSE,
              mvPlot=FALSE)

d_prop05 <- addItem(d_prop05, item = fit025, itemName="ReplicateGroup_fit_eBayes", itemType= "fit")
contrastList  <- list(CDAHFD_vs_Control = "ReplicateGroupMCD.HFD_Veh - ReplicateGroupControl_Veh",
                      BMS963272_3_vs_CDAHFD = "ReplicateGroupMCD.HFD_BMS963272_Hi - ReplicateGroupMCD.HFD_Veh"
)
d_prop05 <- runContrasts(d_prop05,
                       designMatrixName="ReplicateGroup",
                       contrastList=contrastList,
                       runTopTreat=TRUE,
                       FoldChangeThreshold = 2,
                       Qvalue=TRUE,
                       IHW=TRUE,
                       runEBayes=TRUE,
                       proportion=0.05,
                       verbose=TRUE)
cat("***** Proportion = 0.05 Results: *****\n")
cat("No FC Threshold\n")
knitr::kable(summarizeSigCounts(getType(d_prop05, "topTable")))
cat("FC Threshold = 2\n")
knitr::kable(summarizeSigCounts(getType(d_prop05, "topTable"), fcThreshold = 2))
cat("treat FC Threshold = 2\n")
knitr::kable(summarizeSigCounts(getType(d_prop05, "topTreat")))


############################
## Rerun with proportion = 0.25
############################

d_prop25 <- runVoom(d2,
                     designMatrixName="ReplicateGroup",
                     runEBayes = FALSE,
                     mvPlot=FALSE)

contrastList  <- list(CDAHFD_vs_Control = "ReplicateGroupMCD.HFD_Veh - ReplicateGroupControl_Veh",
                      BMS963272_3_vs_CDAHFD = "ReplicateGroupMCD.HFD_BMS963272_Hi - ReplicateGroupMCD.HFD_Veh"
)
d_prop25 <- runContrasts(d_prop25,
                          designMatrixName="ReplicateGroup",
                          contrastList=contrastList,
                          runTopTreat=TRUE,
                          FoldChangeThreshold = 2,
                          Qvalue=TRUE,
                          IHW=TRUE,
                         runEBayes=TRUE,
                         proportion=0.25,
                         verbose=TRUE)

cat("Proportion = 0.25 Results:\n")
cat("No FC Threshold\n")
knitr::kable(summarizeSigCounts(getType(d_prop25, "topTable")))
cat("FC Threshold = 2\n")
knitr::kable(summarizeSigCounts(getType(d_prop25, "topTable"), fcThreshold = 2))
cat("treat FC Threshold = 2\n")
knitr::kable(summarizeSigCounts(getType(d_prop25, "topTreat")))
