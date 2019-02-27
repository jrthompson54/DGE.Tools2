#### Setup ###

library(magrittr)
library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(ArrayServer)
library(zFPKM)

setwd("~/P01380_BuildDGEobj/IBD_FLA")
inputPath <- "./input"
outputPath <- "./output"

pid <- c("P-20170717-0001", "P-20170807-0002")
projectName <- "RNA-Seq_Analysis_of_FL_IBD_Human_Biopsy_P-20170717-0001" #as it appears in Omicsoft and the Regfile
rdsname <- "RNA-Seq_Analysis_of_FL_IBD_Human_Biopsy_P-20170717-0001"
regfile <- file.path(inputPath, str_c(projectName, ".txt"))  #omicsoft registration file; .txt or .gz file
designMatrixName <- "Category2"   #appropriate name for this model
formula <- '~ 0 + Category2'

#Use Design$TRDSample.x for this dataset
dupcorBlock <- NULL    #define duplicates for dupliceCorrelation method; set to NULL to disable

#build contrastList from design matrix column names
contrastList <- list(
  CD_LesionVsNonlesion_Colon =  "CD_Lesion_Colon - CD_Nonlesion_Colon",
  UC_LesionVsNonlesion_Colon = "UC_Lesion_Colon - UC_Nonlesion_Colon",
  CD_LesionVsControl_Colon = "CD_Lesion_Colon - Control_Control_Colon",  
  CD_NonlesionVsControl_Colon = "CD_Nonlesion_Colon - Control_Control_Colon",
  UC_LesionVsControl_Colon = "UC_Lesion_Colon - Control_Control_Colon",  
  UC_NonlesionVsControl_Colon = "UC_Nonlesion_Colon - Control_Control_Colon",
  
  CD_LesionVsNonlesion_Ileum = "CD_Lesion_Ileum - CD_Nonlesion_Ileum",
  UC_LesionVsNonlesion_Ileum = "UC_Lesion_Ileum - UC_Nonlesion_Ileum",
  CD_LesionVsControl_Ileum = "CD_Lesion_Ileum - Control_Control_Ileum",
  CD_NonlesionVsControl_Ileum = "CD_Nonlesion_Ileum - Control_Control_Ileum",
  UC_LesionVsControl_Ileum = "UC_Lesion_Ileum - Control_Control_Ileum",
  UC_NonlesionVsControl_Ileum = "UC_Nonlesion_Ileum - Control_Control_Ileum"
)

#### End of Setup ###
#
# Aside from the Setup block; see lines following a string of ##### for lines
# further down that may need to be modified depending on the specific project.
# line 59, 64, 73
#
#
RDSname <- file.path(outputPath, str_c(rdsname, ".RDS"))


#     #  getFilez grabs Omicsoft data from an S3 bucket.
#     getFilez <- function(projectID, projectName, outputDir="./", 
#                          files=c("RNA-Seq.Design.txt", 
#                                  "RNA-Seq.Count.Table.txt",
#                                  "RNA-Seq.Count.Annotation.txt",
#                                  "RNA-Seq.QCMetrics.Table.txt")) {
#       #Omicsoft data is organized in the S3 bucket by projectID/projectName/ExportedViewAndTables/filename
#       if (missing(projectName)){
#         print("Available Projects:\n")
#         print(getOmicsoftProjectList(projectid=projectID))
#         print("\n")
#         stop()
#       }
#       #save current path
#       cpath <- getwd()
#       setwd(outputDir)
#       getArrayServerFile(projectid=projectID,
#                          omicsoftproject=projectName,
#                          files=files)
#       system("gzip *.txt") #assumes a gzip command is installed and pathed (e.g. cygwin)
#       setwd(cpath)
#     }
# 
# #get Omicsoft data from cloud
# aws.signature::use_credentials()
# ########################################  may or may not be a .gz file
# if (!file.exists(file.path(inputPath, "RNA-Seq.Count.Table.txt.gz"))){
#     getFilez(pid, projectName=projectName,
#          outputDir=inputPath,
#          files=c("RNA-Seq.Design.txt/RNA-Seq.Design.txt", 
#                  "RNA-Seq.Count.Table.txt",
#                  "RNA-Seq.Count.Annotation.txt/RNA-Seq.Count.Annotation.txt",
#                  "RNA-Seq.QCMetrics.Table.txt"))
# }
# ######################################### set gz values as needed
# dgeObj <- OmicsoftToDgeObj(path=inputPath, gz=TRUE)

#####  For Xpress Data  #######
# uncomment lines with single # to use this block to retrive data from Xpress
# You'll need the xpress project number which is the tail of the URL for that project  e.g:
# http://xpress.pri.bms.com/CGI/project_summary.cgi?project=19958
# The Xpress project ID is 19958
XpressID = 20245
library(Xpress2R)
# First return the available "levels"  (typically genes or isoforms but sometimes multiple genes with different gene models)
# > dgeObj <- Xpress2DGEO(XpressID)
# Message: mm10MTERCC-ensembl80-genes, mm10MTERCC-ensembl80-isoforms
# Now we know what to ask for
if (!file.exists(file.path(outputPath, "Xpress20245.RDS"))){
  dgeObj <- Xpress2DGEO(XpressID, level= "GRCh37ERCC-ensembl75-genes") #might take 2-3 min
  saveRDS (dgeObj, file.path(outputPath, "Xpress20245.RDS"))
} else {
  dgeObj <- readRDS(file.path(outputPath, "Xpress20245.RDS" ))
}
# 
#####  End Xpress Data Block #####
# dgeObj <- readRDS(file.path(outputPath, "IBD_FLA.RDS" ))

#add project metadata
dgeObj <-  annotateDGEobj(dgeObj, regfile=regfile)

#import new design table
newDesign <- read.delim("./input/CombinedDOE_ProteomicsNAseq_257_noCrossValSamples.txt",
                                 stringsAsFactors=FALSE)
rownames(newDesign) <- str_replace(newDesign$smps, "GRCh37ERCC_ensembl75_genes_", "")
#Trim the DGEobj to samples in newDesign
idx <- colnames(dgeObj) %in% rownames(newDesign)
dgeObj <- dgeObj[,idx]
#correct order in newDesign
all(rownames(newDesign) %in% colnames(dgeObj))
newDesign <- newDesign[colnames(dgeObj),]
all(rownames(newDesign) == colnames(dgeObj))


#create ReplicateGroup column
newDesign$ReplicateGroup <- newDesign$Category2


dgeObj <- addItem(dgeObj, newDesign, itemName = "design", itemType="design", overwrite=TRUE)
dim(dgeObj) #57905 rows

#Filter out genes with zero efflength (only needed for data from Xpress)
el <- getItem(dgeObj, "effectiveLength")
idx <- el == 0
el[idx] <- NA
idx <- complete.cases(el)
dgeObj <- dgeObj[idx,]
dim(dgeObj) #50534 rows

#Low Intensity Filter
fracThreshold <- 0.5

dgeObj <- lowIntFilter(dgeObj, zfpkmThreshold = -3, countThreshold = 10, sampleFraction=fracThreshold)

#low expression filter
# #keep zfpkm >= -3 in fracThreshold of samples
# fpkm <- convertCounts(dgeObj$counts, unit="fpkm", log=FALSE, normalize = "tmm", 
#                       geneLength=dgeObj$geneData$ExonLength) %>% as.data.frame
# zfpkm <- zFPKM(fpkm)
# idxzfpkm <- zfpkm >= -3.0
# frac <- rowSums(idxzfpkm) / ncol(idxzfpkm)
# idx <- frac >= fracThreshold
# dgeObj <- dgeObj[idx,]
# dim(dgeObj) #25218 Rows
# 
# #overlay a mincount filter
# counts <- getItem(dgeObj, "counts")
# # genelength <-getItem(dgeObj, "geneData")$ExonLength
# idxmin <- counts >= 10
# frac <- rowSums(idxmin)/ncol(idxmin)
# idx <- frac >= fracThreshold
# dgeObj <- dgeObj[idx,] #19267
# 
# #Limit to protein coding
# gd <- dgeObj$geneData
# idx <- dgeObj$geneData$Biotype == "protein_coding"
# idx[is.na(idx)] <- FALSE  #mt genes have NA bioptype
# dgeObj <- dgeObj[idx,]
# dim(dgeObj)

# #tempsave
# saveRDS(dgeObj, file.path(outputPath, "IBD_FLA_filtered.RDS"))

# Build a design matrix and store in DGEobj
design <- getItem(dgeObj, "design")


#Make model terms into factors (set 1st level to desired baseline)
##########################################
##########################################

# design$Category2 <- factor(design$ReplicateGroup, levels=c("naive_veh_D1", "Sham_veh_D9", "UUO_veh_D9", "UUO_EXT001024_D9"))

#set up and save the design matrix
designMatrix <- model.matrix (as.formula(formula), design)

colnames(designMatrix) <- str_remove(colnames(designMatrix), "Category2")
#can't have colons in colnames
colnames(designMatrix) <- make.names(colnames(designMatrix))


#capture the formula as an attribute
designMatrix <- setAttributes(designMatrix, list(formula=formula))
#save the modified designMatrix
dgeObj <- addItem(dgeObj, item=designMatrix, 
                  itemName=designMatrixName, 
                  itemType="designMatrix",
                  parent="design", overwrite=TRUE)


### run DGE calculations
#normalize
dgeObj <- runEdgeRNorm(dgeObj, plotFile=file.path(outputPath, str_c("TMM_Norm.Factors.PNG")), normMethod = "TMM")  #TMM is the default

dupcorBlock <- dgeObj$design$TRDSample.x

print("Running Voom")

#voom/lmfit
if (is.null(dupcorBlock)){
  dgeObj <- runVoom(dgeObj, designMatrixName,
                    qualityWeights = TRUE,
                    mvPlot = FALSE)
} else {
  dgeObj <- runVoom(dgeObj, designMatrixName,
                    qualityWeights = TRUE,
                    dupcorBlock = dupcorBlock)
}

print("running contrasts")

#run contrasts
dgeObj <- runContrasts(dgeObj, 
                       designMatrixName=designMatrixName, 
                       contrastList=contrastList, 
                       runTopTreat=FALSE,
                       Qvalue=TRUE,
                       IHW=TRUE)

saveRDS(dgeObj, RDSname)

attr(dgeObj, "source") <- "IBD_FLA_DGEpipeline.R"
attr(dgeObj, "repo") <- "https://biogit.pri.bms.com/search?utf8=%E2%9C%93&q=P01380_BuildDGEobj&type="
saveRDS(dgeObj, RDSname)

knitr::kable(inventory(dgeObj))

knitr::kable(summarizeSigCounts(getType(dgeObj, "topTable"), fcThreshold = 1.5), caption="FLIBD without DupCor")

JRTutil::checkDGEobj(dgeObj)
              