### Function Plot Dispersion ###
#' Function  plotDispersion
#'
#' Creates an edgeR dispersion plot for RNA-Seq QC purposes.  Takes a SLOA object
#' for input after runEdgeRNorm has been executed.  The Biological coefficient
#' of variation (BCV) is displayed.  BCV is the sqrt of dispersion.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, dispersion, plot, QC
#'
#' @param MySLOA A SLOA object from DGE:Tools runEdgeRNorm with runDisp=T (the
#'   default).
#' @param symbolSize (Default=1)
#' @param symbolShape see
#'   http://www.cookbook-r.com/Graphs/Shapes_and_line_types/ (Default = 1)
#' @param symbolColor Default = "darkblue"
#' @param symbolFill  Default = "darkblue"
#' @param symbolAlpha Transparency for the points. Value from 0 to 1. Smaller
#'   indicate more transparency (Default = 0.3)
#' @param linefitSize (Default = 1)
#' @param linefitColor (Default = "yellow")
#' @param lineFit (Default = NULL)  Defaults to NULL because enabling this
#'   overloads my poor 32bit PC. Any type supported by geom_smooth (loess is
#'   recommended).
#' @param rugColor (Default = NULL)  Set to Null to disable rug plots.
#'   Note: printing plots with rug plots can be slow (several minutes on my
#'   poor 32bit pc). So I left this diabled by default.  Assign a valid
#'   color to enable
#' @param rugAlpha Transparency for the rug plots (Default = 0.02)
#'
#' @return a ggplot object
#'
#' @examples
#'    MyGgplot <- plotDisp (MyData=MySLOA)
#'    MyGgplot <- plotDisp (MyData=Mydgelist, DesignMatrix=MyDesignMatrix)
#'
#' @import edgeR
#'
#' @export
plotDisp <- function (MySLOA,
                            symbolSize = 1,
                            symbolShape = 1,
                            symbolColor = "darkblue",
                            symbolFill = "darkblue",
                            symbolAlpha = 0.3,
                            linefitSize = 1,
                            linefitColor = "yellow",
                            lineFit = NULL,
                            lineType = 1,
                            lineAlpha = 1,
                            rugColor = NULL,
                            rugAlpha = 0.02) {
# browser()
  #dispersion: expect 0.4 coef of var for human, 0.1 for inbreed mouse, near 0 for technical reps.
  #  source("~/R/lib/SubsettableListOfArrays.R")

  #Check for and source SubsettableListOfArrays.R  Will save people from
  #remembering to source it manually if they put it where I suggested.
  if (file.exists("~/R/lib/SubsettableListOfArrays.R")){
    source("~/R/lib/SubsettableListOfArrays.R")
  }

  if (!exists("SubsettableListOfArrays")){
    error("You must source SubsettableListOfArrays.R first.  http://bioinformatics.bms.com/active/biohtml/thompj27/DGE.Tools/SubsettableListOfArrays.R.")
  }

  if (!class(MySLOA)[[1]] == "SLOA"){
      #get the required data from the SLOA object
      error("MySLOA must be class SLOA")
  }

  plotdata <- data.frame(AveLogCPM=MySLOA$DGElist$AveLogCPM, BCV=sqrt(MySLOA$DGElist$tagwise.dispersion))

  dgelist <- MySLOA$DGElist
  MyDispPlot <- ggplot(plotdata, aes(x=AveLogCPM, y=BCV)) +
    geom_point (size=symbolSize, shape=symbolShape, fill=symbolFill,
                color=symbolColor, alpha=symbolAlpha)

  if (!is.null(lineFit)){
   MyDispPlot <- MyDispPlot +
      geom_smooth(method=lineFit, size=linefitSize, color=linefitColor,
                linetype=lineType, alpha=lineAlpha)
  }

  if (!is.null(rugColor)){
    MyDispPlot = MyDispPlot + geom_rug(data=plotdata, inherit.aes=FALSE,
                                         color = rugColor,
                                         alpha=rugAlpha,
                                         show.legend=FALSE,
                                         aes(x=AveLogCPM, y=BCV))
  }

  MyDispPlot <- MyDispPlot +
    ggtitle ("EdgeR Dispersion Plot") +
    greyTheme()

  return(MyDispPlot)
}
