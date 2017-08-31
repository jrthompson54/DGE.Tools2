library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(Xpress2R)
library(magrittr)
library(tidyverse)
library(stringr)
#Get CCL4 data from Omicsoft
do <- readRDS("~\\R\\contrastDB\\DGEobj\\DGEobj_APJ_LAD-heart_totalRNA(P-20170103-0001).RDS")
do <- annotateDGEobj(do, "~\\R\\contrastDB\\DGEobj\\P-20170103-0001_totalRNA_Registration.txt")
attr(do, "source") <- "Omicsoft"
attr(do, "PID") 
names(attributes(do))
dim(do)



#Get XpressData
#http://xpress.pri.bms.com/CGI/project_summary.cgi?project=20140
#PID P-20170103-0001
if (file.exists("X20140.RDS")) {
    dx <- readRDS("X20140.RDS")
} else {
    dx <- Xpress2DGEO(20135)
    saveRDS (dx, "X20140.RDS")
}

#Get Xpress TPM data
library(rXpress)
projectNumber <- 20140
variables <- rXpress::getVariables(projectNumber)
variables$projectInfo$ATTRIBUTES
attributes <- getAttributes(projectNumber)
xpressData <- rXpress::getXpressData(projectNumber, attribute=201)
varSetNames <- rXpress::getVarSetNames(projectNumber)
Level <- varSetNames[1]
tpm.xpress = xpressData$varSets[[Level]]$y
###Filter out non-Ensembl records 
idx <- str_detect(rownames(tpm.xpress), "^ENS")  #ENSG or ENST ensembl records
tpm.xpress <- tpm.xpress[idx,]
#Remove genome reference from colnames
colnames(tpm.xpress) <- str_sub(colnames(tpm.xpress), (nchar(Level)+2))


#print length==1 attributes in a table
x <- attributes(do)
idx <- lapply(x, length) == 1
y <- stack(x[idx])
colnames(y) <- c("Value", "Attribute")
y <- select(y, Attribute, Value)
y

#print length==1 attributes in a table
x <- attributes(dx)
idx <- lapply(x, length) == 1
y <- stack(x[idx])
colnames(y) <- c("Value", "Attribute")
y <- select(y, Attribute, Value)
y


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
#looks consistent.  


#now try with Xpress
tpm.dx <- tpm(dx)
dim(tpm.dx)
dim(counts.dx)

#dataset was floored!  Crudely reverse it.
idx <- counts.dx == 0.01
counts.dx[idx] <- 0
tpm.dx[idx] <- 0
inspect(counts.dx)
inspect(tpm.dx)
min(counts.dx)
#counts are still being floored at 0.01!
# idx <- counts.dx == 0.01
# counts.dx[idx] <- 0
inspect(counts.dx)

#note 48% zeros in Omicsoft data; 65% zeros in Xpress data


# Merge tpm from xpress and xpress via tpm function for the first sample
test <- cbind(tpm.dx[,1], tpm.xpress[,1]) %>% as.data.frame
