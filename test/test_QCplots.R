library(tidyverse)
library(DGEobj)
library(DGE.Tools2)

s3mount <- "Y:"
s3path <- "/OmicsoftHome/output/P-20180326-0001/TempleUniv_HeartFailure2017_P-20180326-0001_R94_24Jan2019/ExportedViewsAndTables"
qcfilename <- "RNA-Seq.QCMetrics.Table.txt"
qcdat <- readr::read_delim(file.path(s3mount, s3path, qcfilename), delim="\t")
colnames(qcdat) <- stringr::str_sub(colnames(qcdat), 1, 16)


metrics <- c("Alignment_MappedRate", "Alignment_PairedRate", "Source_rRNA", "Strand_Read1AntiSenseRate",
             "Strand_ReadPairAntiSenseRate", "Profile_ExonRate", "Profile_InterGene_FPK")
m <- "Alignment_MappedRate"
plots <- QCplots(qcdat, metricNames=metrics, plotType="bar", xAngle=90)
plots_NoWin <- QCplots(qcdat, metricNames=metrics, plotType="bar", xAngle=90, winsorize=FALSE)

plots <- QCplots(qcdat, metricNames=m, plotType="bar", xAngle=45)

x <- ggplot(qcdata, aes(x=Sample, y=Alignment_MappedRate))+
    geom_bar(stat="identity", fill="dodgerblue3")

yvar <- "Alignment_MappedRate"
x <- ggplot(qcdata, aes(x=Sample, y=!!!yvar))+
  geom_bar(stat="identity", fill="dodgerblue3")
