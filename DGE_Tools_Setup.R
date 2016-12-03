#DGE_Tools Setup
#Restart R and don't load any libraries before running this script (loaded libraries can't be updated)


Install.BiocLite.diff <- function (dependencies) {
#I've noticed that using install_git seems to install dependencies from CRAN
#but not from Bioconductor.
#
#Pass a vector of dependencies to this funtion
#The function queries for a list of user installed packages, performs
#a diff and attempts to install the missing packages with biocLite.
#
  source("https://bioconductor.org/biocLite.R")

  #get user installed packages
  #http://www.r-bloggers.com/list-of-user-installed-r-packages-and-their-versions/
  ip <- as.data.frame(installed.packages()[,c(1,3:4)])
  rownames(ip) <- NULL
  ip <- ip[is.na(ip$Priority),1:2,drop=FALSE]
  uinstalled = ip$Package

  missing = setdiff(dependencies , uinstalled)

  biocLite(missing, suppressUpdates=TRUE)
}

dependencies = c(
  "AnnotationDbi",
  "dplyr",
  "edgeR",
  "gdata",
  "GenomicRanges",
  "ggplot2",
  "grDevices",
  "grid",
  "gridExtra",
  "IRanges",
  "limma",
  "locfit",
  "magrittr",
  "openxlsx",
  "org.Hs.eg.db",
  "org.Mm.eg.db",
  "org.Rn.eg.db",
  "pheatmap",
  "qvalue",
  "qvalue",
  "RColorBrewer",
  "reshape2",
  "S4Vectors",
  "stringi",
  "stringr",
  "SummarizedExperiment",
  "tidyr",
  "tools"
)


Install.BiocLite.diff(dependencies)

#Install DGE.Tools and zFPKM
#presumes you've already copied DGE.Tools_0.0.0.9000.tar.gz and zFPKM_0.0.0.9000.tar.gz to ~/R/lib
#install.packages("~/R/lib/zFPKM_0.0.0.9000.tar.gz", repos=NULL, type="source")
#install.packages("~/R/lib/DGE.Tools_0.0.0.9011.tar.gz", repos=NULL, type="source")

install_git("http://biogit.pri.bms.com/thompj27/zFPKM")
install_git("http://biogit.pri.bms.com/thompj27/DGE.Tools")

print("Don't forget to download SubsettableListOfArrays.R too")

#fail
#install_github ("thompj27/R-Packages/DGE.Tools/prod", host="http://biogit.pri.bms.com")
#install_git("http://biogit.pri.bms.com/thompj27/R-Packages/tree/master/DGE.Tools/prod")

#I think the un/pw dialog is killing the download this way (invalid octal digit)
# install.packages(
#  "http://bioinformatics.bms.com/active/biohtml/thompj27/DGE.Tools/DGE.Tools_0.9014.tar.gz",
#  repos = NULL, type = "source"
# )

#  http://bioinformatics.bms.com/active/biohtml/thompj27/DGE.Tools/DGE.Tools_0.9014.tar.gz
#

# wdsave = getwd()
# setwd("~/R/lib")
# download.file(
#   "http://bioinformatics.bms.com/active/biohtml/thompj27/DGE.Tools/DGE.Tools_0.9014.tar.gz",
#   "DGE.Tools_0.9014.tar.gz"
# )
# install.packages("DGE.Tools_0.9014.tar.gz", repos = NULL, type = "source")
# setwd(wdsave)
