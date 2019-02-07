library(tidyverse)
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)

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

plots <- QCplots(qcdat, metricNames=metrics, plotType="histogram", xAngle=0)


#test JRTutil winsorize functions
qcdata <- qcdat %>% column_to_rownames(var="Metric") %>% t() %>% as.data.frame

missingDataHeatmap(qcdata)

mean(qcdata[,3,drop=T], na.rm=T)
winsorize_mean(qcdata[,3,drop=T])
median(qcdata[,3,drop=T])

sd(qcdata[,3,drop=T])
winsorize_sd(qcdata[,3,drop=T])

var(qcdata[,3,drop=T])
winsorize_var(qcdata[,3,drop=T])



#QC heatmap scaling each metric.
#380 cells =NA
idx <- complete.cases(qcdat)
qcdat <- qcdat[idx,]
qcdata <- qcdat %>% column_to_rownames(var="Metric") %>% t() %>% as.data.frame
pheatmap::pheatmap(qcdata, scale = "column", breaks=NA)
test <- qcdata[!idx,]

missingDataHeatmap(qcdata)
inspect(qcdata)
