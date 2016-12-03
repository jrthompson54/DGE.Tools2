#JRT 26Sep2015

#' @include DGE_Tools.R

#use Build_RSE_Object to run Build_RSE and save the result to a specified file.
#Can be run from the R command line
#Designed to by run from the DOS commandline by calling Rscript Run_Build_RSE.R
#
#Defaults to running GeneData level with ./ as the input dir, ./RData as the
#output folder, zFPKM.PNG as the plotfile and RSE.RDS as the output file,
#facet titles (on the zFPKM plot) turned off.
#

#suppressMessages({
#
#  library(magrittr)
#
#  #Meaning of ~ is different at the DOS commandline and within RStudio
#  #In RStudio ~ = c:/Users/jrt/Documents
#  #In DOS ~ = c:/Users/jrt

#  if (file.exists("~/Documents/R/lib/DGE_Tools.R")){
#    source("~/Documents/R/lib/DGE_Tools.R", chdir = TRUE)
#  } else if (file.exists("~/R/lib/DGE_Tools.R")) {
#    try(source("~/R/lib/DGE_Tools.r", chdir = TRUE))
#  } else if (file.exists("/home/thompj27/R/lib/DGE_Tools.R")) {
#    source("/home/thompj27/R/lib/DGE_Tools.R", chdir = TRUE)
#  } else {
#    stop(sprintf ("Can't file DGE_Tools.R!  Halting."))
#  }
#  library(stringr)
#})

trim_trailing_pathsep = function(txt) {
  #removes trailing slash or backslash or any string specified by platformsep
  pathsep.regex = sprintf("[%s\\\\/]+$", .Platform$file.sep)
  str_replace(txt, pathsep.regex, "")
}

### Function BuildAndSaveRSE ###
#' Function  BuildAndSaveRSE
#'
#' Build a RangedSummarizedExperiment (RSE) from Omicsoft Alignment and count data.  The
#' RSE is saved as a .RDS file.  Users provides
#' separate tab-delimited text files for each row annotation, column anotation
#' and "assay" (matrices of RxC e.g. counts, FPKM). GeneData is a predefined data
#' structure (list of lists) that specifies the datafiles to import.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID, SummarizedExperiment, Omicsoft
#'
#' @param InputFolder defaults to current folder
#' @param OutputFolder defaults to "./RData"
#' @param RSEfilename Filename for the RSE.RDS output file
#' @param PlotFilename Filename for the zFPKM PNG file.
#' @param Species One of Human, Mouse or Rat.
#' @param Level One of "GeneData" or "TranscriptData" (defaults to "GeneData")
#' @param FacetTitles  Name for image of sample distribution fits
#'
#' @return Datafram with Entez GeneID column added
#'
#' @examples
#' MyDataframe = Build_RSE (GeneData, Species = "Human", OutputPath = "./RData", PlotFile = "Gene.zFPKM.PNG")
#'
#' @import zFPKM dplyr SummarizedExperiment GenomicRanges
#'
#' @export
BuildAndSaveRSE <- function(InputFolder = ".", OutputFolder = "./RData", RSEfilename = "RSE.RDS",
  Species = "", PlotFilename = "zFPKM.PNG", Level="GeneData", FacetTitles = FALSE) {
# JRT 3Sep2015
# BuildAndSaveRSE imports various RNA-Seq output files from the Omicsoft pipeline and
# creates a RangedSummarizedExperiment that is saved in an .RData file.
#
# Primarily intended to be run by the command line script Run_Build_RSE.R
#
# Input parameters are all optional.
# Build_RSE_Object() will create the file ./RData/RSE.RData
# RSE.RData will contain a single variable (named RSE) which is a RangedSummarizedExperiment
# and will contain Gene level data.
#
# Use the OutputFolder and RSEfilename parameters to change the output folder and the output
# filename respectively.  The output filename should always have a .RData extension.
#
# Use Level="TranscriptData" to operate on transcript level data instead of gene level data.
#
# See the GeneData and Transcript Data lists defined in DGE_Tools.r to see the text files
# that must be present for this to operate.
#
# Usage:
# Copy this script into the folder containing your RNA-Seq data files
# Then execute this script from a terminal with the command:
#     RScript Build_RSE_Object.R

InputFolder = trim_trailing_pathsep(InputFolder)
OutputFolder = trim_trailing_pathsep(OutputFolder)

  if (!file.exists(InputFolder)) {
    stop (sprintf ("Input Directory \"%s\" does not exist. Halted!", InputFolder))
  } else {
    setwd(InputFolder)
  }

  #Set the level
  if (tolower(Level) == "transcriptdata") {
    RSE = Build_RSE(TranscriptData, PlotFile=PlotFilename, FacetTitles=FacetTitles,
                    Species=Species, OutputPath=OutputFolder)
  } else if (tolower(Level) == "genedata") {
    RSE = Build_RSE(GeneData, PlotFile=PlotFilename, FacetTitles=FacetTitles,
                    Species=Species, OutputPath=OutputFolder)
  } else {
    stop ("Error: Level must be \"TranscriptData\" or \"GeneData\"")
  }

  #create the output folder, if needed
  if (!file.exists(OutputFolder)){
    dir.create(OutputFolder)
  }

  saveRDS(RSE, file=file.path(OutputFolder, RSEfilename))
}

