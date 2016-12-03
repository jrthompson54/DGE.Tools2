### Function Plot Dispersion ###
#' Function  plotDispersion
#'
#' Warning: Deprecated. Kept for now for backward compatibility but will be
#' removed in a later version. Replaced by DGE.Tools::plotDisp. See ?plotDisp
#' for details on the replacement function.
#'
#' Runs an edgeR dispersion plot for RNA-Seq QC purposes.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, dispersion, plot
#'
#' @param MyCounts A martrix of count data with geneID as rownames.
#' @param MyFormula a text string specifying the model formula
#' @param MyDesign  Experiment design dataframe with columns referenced by the forumula
#' @param showPlot  Set TRUE to print the plot
#' @param savePlot Set TRUE to save plot in a .PNG file
#'
#' @return nothing
#'
#' @examples
#'    plotDispersion (MyCounts, MyDesignMatrix, showPlot=TRUE, savePlot=TRUE)
#'
#' @import edgeR
#'
#' @export
plotDispersion <- function (MyCounts, MyFormula=NULL, MyDesign=NULL, showPlot=TRUE, savePlot=TRUE) {
  # Input:
  #   MyCounts is a Counts matrix
  #   MyDesignMatrix is a design matrix created by model.matrix function.
  #   Optionally, showPlot and savePlot can be used to suppress one or the other of the plots
  #
  # Output:  Draws a plot (dependent on showPlot) and optionally saves the plot
  #          to a file (dependent on savePlot)


  FVersion = "plotDispersion : 26Mar2016"

  #dispersion: expect 0.4 coef of var for human, 0.1 for inbreed mouse, near 0 for technical reps.
  #  source("~/R/lib/SubsettableListOfArrays.R")

  #Build design matrix from formula and design data.frame
  #
  warning("Warning: plotDispersion is Deprecated. Replaced by DGE.Tools::plotDisp. See ?plotDisp for details")

  MyDesignMatrix = model.matrix(as.formula(MyFormula), MyDesign)

  dgeList = MyCounts %>%
    as.matrix %>%
    DGEList %>%
    calcNormFactors %>%
    estimateDisp (design=MyDesignMatrix, robust=TRUE)

  if (showPlot) {
    invisible(plotBCV(dgeList, ylab="Biol. Coeff. of Variation", cex.axis=1.5, cex.lab=1.5))
  }
  if (savePlot) {
    png(file="Dispersion.png",width=8,height=6, units = 'in', res = 300)
    plotBCV(dgeList, ylab="Biol. Coeff. of Variation", cex.axis=1.5, cex.lab=1.5)
    invisible ( dev.off() )
  }

  return(dgeList)
}
