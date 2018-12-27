### Function logRatioPlot ###
#' Function logRatioPlot
#'
#' Intended to plot a set of contrast results, one plot for each gene of
#' interest. Input is a tidy datafile constructed from topTable output and
#' requires logFC, CI.L and CI.R columns as well as a gene identifier of choice.
#' Outputs a ggplot object facetted by the facetColname or a list of individual
#' ggplots, one for each facetColname value (typically gene).
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords barplot lineplot gplot2 logratio confidence intervals contrasts
#'
#' @param data A tidy dataframe of data to plot (see ?tidyContrasts).
#' @param facetColname Define the column name to separate plots (e.g. GeneID).
#' @param xColname Define the column name to group boxplots by (e.g. Constrast).
#' @param yColame Define the column of values for plotting (Default = "logFC").
#' @param CI.R_colname Define name of the CI high value (Default = "CI.R")
#' @param CI.L_colname Define name of the CI low value (Default =  "CI.L")
#' @param xOrder Define the order for the groups in each plot.  Should
#'   contain values in unique(data$group) listed in the order that you want the
#'   groups to appear in the plot. (optional; default = unique(data[xColname]))
#' @param plotType one of "bar" or "point" (Default = bar")
#' @param addLine Adds a line if point layer chosen (Default = FALSE)
#' @param refLine Adds a horizontal line at y=0 (Default=TRUE)
#' @param refLineColor Color for the reference line (Default="red")
#' @param xlab X axis label (defaults to xColname)
#' @param ylab Y axis label (defaults to yColname)
#' @param title Plot title (optional)
#' @param barColor Color for the bar outline (default = "dodgerblue4")
#' @param barFill Color for the bar area (default = "dodgerblue3")
#' @param barSize set the bar size (thickness of each bar perimeter; default = 0.1)
#' @param barWidth set the bar width (Default = 0.8)
#' @param barAlpha Transparency for the bar layer (Default = 1)
#' @param pointColor Color for the point layer (Default = "grey30")
#' @param pointFill Fill color for the point layer (Default = "dodgerblue4")
#' @param pointShape Shape for the point layer (Default = 21; fillable circle)
#' @param pointAlpha Transparency for the box layer (Default = 1)
#' @param pointSize Size of the points (Default = 4)
#' @param lineLayer Add a fitted line layer (Default = FALSE)
#' @param lineColor Color of the line fit (Default = "dodgerblue4")
#' @param lineSize Size of the line fit (Default = 1)
#' @param lineFit Type of fit to use.  One of c("auto", "lm", "glm", "gam",
#'   "loess"). (Default = "loess")
#' @param lineType One of c("solid", "dashed", "dotted", "dotdash", "longdash",
#'   "twodash"). (Default = "solid")
#' @param baseFontSize The smallest size font in the figure in points. (Default =
#'   12)
#' @param themeStyle "bw" or "grey" which correspond to theme_bw or theme_grey
#'   respectively. (Default = "bw")
#' @param facet Specifies whether to facet (TRUE) or print individual plots
#'   (FALSE)  (Default = TRUE)
#' @param facetCol Explicitly set the number of Rows for the facet plot. Default
#'   behavior will automatically set the columns. (Default = ceiling(sqrt(length(unique(data[facetCol])))))
#' @param xAngle Angle to set the sample labels on the Xaxis (Default =  45; Range = 0-90)
#' @param scales Specify same scales or independent scales for each subplot (Default = "free_y";
#'   Allowed values: "fixed", "free_x", "free_y", "free")
#' @param debug Turn on debug mode (default = FALSE)
#'
#' @return ggplot If Facet=TRUE (default) returns a facetted ggplot object. If
#'   facet=FALSE, returns a list of ggplot objects indexed
#'   by observation (gene) names.
#'
#' @examples
#'
#'   #DGEobj example
#'   contrastList <- getType(dgeObj, "topTable")
#'   #Put contrasts in tidy format keeping logFC, and confidence limits data
#'   tidyDat <-tidyContrasts(dgeObj, rownameColumn="EnsgID", includeColumns = c("logFC", "CI.R", "CI.L")
#'   #add gene symbols from geneData
#'   ens2genesym <- dgeObj$geneData %>%
#'                  rownames_to_column(var="EnsgID") %>%
#'                  select(EnsgID, GeneSymbol=GeneName)
#'  tidyDat <-  left_join(tidyDat, ens2genesym)
#'  #simple barplot
#'  logRatioPlot(tcDat,
#'               facetColname = "GeneSymbol",
#'               xColname = "Contrast")
#'
#'  #lineplot with some options
#'  logRatioPlot(tcDat, plotType="point",
#'                      facetColname = "GeneSymbol",
#'                      xColname = "Contrast",
#'                      facetCol=4,
#'                      scales="fixed",
#'                      facet=TRUE,
#'                      title = "Test",
#'                      pointSize=4,
#'                      lineLayer=TRUE,
#'                      lineSize=0.1,
#'                      xAngle=60)
#'
#' @import ggplot2 magrittr
#' @importFrom assertthat assert_that
#' @importFrom stringr str_c
#'
#' @export
logRatioPlot <- function(data,
                      facetColname,
                      xColname,
                      yColname = "logFC",
                      CI.R_colname = "CI.R",
                      CI.L_colname = "CI.L",
                      xOrder=unique(as.character(data[xColname,, drop=TRUE])),
                      plotType = "bar",
                      refLine = TRUE,
                      refLineColor = "red",
                      xlab=xColname, ylab=yColname, title,
                      barColor = "dodgerblue4",
                      barFill = "dodgerblue3",
                      barSize = 0.1,
                      barAlpha = 1,
                      barWidth = 0.9,
                      pointColor = "grey30",
                      pointFill = "dodgerblue4",
                      pointShape = 21, #fillable circle
                      pointAlpha = 1,
                      pointSize = 2,
                      lineLayer = FALSE,
                      lineColor = "dodgerblue4",
                      lineSize = 1,
                      lineType = "solid",
                      lineFit = "loess",
                      baseFontSize = 12,
                      themeStyle = "grey",
                      facet = TRUE,
                      facetCol,
                      xAngle = 45,
                      scales = "free_y",
                      debug = FALSE
                      )
{

      .addGeoms <- function(MyPlot)
        #note uses global values except for MyPlot
      {

        if (tolower(plotType)=="bar"){
          MyPlot = MyPlot + geom_bar(stat="identity",
                                     alpha=barAlpha,
                                     color=barColor,
                                     fill=barFill,
                                     size=barSize,
                                     width=barWidth
          )
        } else if (tolower(plotType)=="point"){
          MyPlot <- MyPlot + geom_point(alpha=pointAlpha,
                                        color=pointColor,
                                        fill=pointFill,
                                        size=pointSize,
                                        shape=pointShape
                                        )
        } else {
          stop ("plotType must be one of c(\"bar\", \"point\")\n")
        }

        #Add error bars if columns present
        if (all(c(CI.L_colname, CI.R_colname) %in% colnames(tcDat) )) {
          MyPlot <- MyPlot + geom_errorbar(aes_string(ymin=CI.L_colname, ymax=CI.R_colname), width=.2)
        } else {
          warning ("Confidence limits columns not found.")
        }

        if (lineLayer == TRUE)
          MyPlot <- MyPlot + geom_smooth(aes_string(group=facetColname),
                                         method=lineFit,
                                         color=lineColor,
                                         size=lineSize,
                                         se=FALSE)

        return(MyPlot)
      }
      ### end .addGeoms

  ### Argument checks

  assertthat::assert_that(!missing(data),
                          class(data)[[1]] == "data.frame",
                          facetColname %in% colnames(data),
                          xColname %in% colnames(data),
                          yColname %in% colnames(data),
                          all(xOrder %in% as.character(data[xColname,,drop=TRUE]))
  )

  if (debug ==TRUE) browser()

### Plot code here
  if (facet) {

    #set facet columns to sqrt of unique observations (rounded up)
    if (missing(facetCol)) {
      numcol <- data[facetCol] %>% unique %>% length %>% sqrt %>% ceiling
    } else {
      numcol = facetCol
    }

    MyPlot <- ggplot2::ggplot (data, aes_string(x=xColname, y=yColname))
    MyPlot <- .addGeoms(MyPlot)
    facetFormula <- stringr::str_c("~", facetColname, sep=" ")
    MyPlot <- MyPlot + ggplot2::facet_wrap(facetFormula, ncol=numcol, scales=scales)

    MyPlot <- MyPlot + ggplot2::xlab(xlab)
    MyPlot <- MyPlot + ggplot2::ylab(ylab)
    if (!missing(title)) MyPlot <- MyPlot + ggplot2::ggtitle(title)
    if (tolower(themeStyle) == "bw" ){
      MyPlot <- MyPlot + theme_bw() + baseTheme(baseFontSize)
    } else {
      MyPlot <- MyPlot + theme_grey() + baseTheme(baseFontSize)
    }

    #rotate xaxis group labels
    if (xAngle > 0){
      MyPlot <- MyPlot + theme(axis.text.x = element_text(angle = xAngle, hjust = 1))
    }

    #Add refLine at 0
    if (refLine == TRUE) MyPlot <- MyPlot + geom_hline(yintercept=0, color=refLineColor, size=0.1)

  } else { #individual plots for each Gene returned in a list

      plotlist <- list()

      for (obs in unique(data[[facetColname]])) {  #for each gene

        dat <- data[data[[facetColname]] == obs, ] #pull data for one gene
        # dat <- dplyr::filter(data, !!facetColname == !!obs)

        aplot <- ggplot(dat, aes(x=xColname, y=yColname)) + #Samples vs Log2CPM
          xlab(xlab) +
          ylab(ylab) +
          ggtitle(obs) +
          theme_grey() + facetTheme(baseFontSize)
        aplot <- .addGeoms(aplot)

        if (!missing(title)) aplot <- aplot + ggplot2::ggtitle (str_c(title, ": ", obs))

        #rotate xaxis group labels
        if (xAngle > 0){
          aplot <- aplot + theme(axis.text.x = element_text(angle = xAngle, hjust = 1))
        }

        #Add refLine at 0
        if (refLine == TRUE) aplot <- aplot + geom_hline(yintercept=0, color=refLineColor, size = 0.1)

        plotlist[[obs]] <- aplot

      }

      MyPlot = plotlist

  }

  return(MyPlot)
}

