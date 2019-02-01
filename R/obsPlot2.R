### Function obsPlot2 ###
#' Function obsPlot2
#'
#' Provides a summary plot (box/violin/points)for each observation (gene), showing data for each
#' experiment group. The plot can optionally include one or more of the
#' following layers: boxplot, violin plot, individual points and/or mean of all
#' points.  The layers are built up in the order listed with user settable
#' transparency, colors etc.  By default, the Boxplot, Point and Mean layers are
#' active. Also, by default, the plots are faceted.  Facet plot can be turned
#' off to return a list of individual ggplot graphics for each gene.
#'
#' Input is a tidy dataframe of intensity data.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords boxplot violinplot ggplot2 logratio
#'
#' @param data A tidy dataframe of intensity data with row and colnames (required)
#' @param plotByCol Define the column name to separate plots (typically a gene ID column) (required)
#' @param groupCol Define the column name to group boxplots by (typically a replicate group column) (required)
#' @param valueCol Define the column of values for plotting (typically a log intensity measure) (required)
#' @param groupOrder Define the order for the groups in each plot.  Should
#'   contain values in unique(data$group) listed in the order that you want the
#'   groups to appear in the plot. (optional; default = unique(data[groupCol]))
#' @param boxLayer Adds a boxplot layer (Default = TRUE)
#' @param violinLayer Adds a violin layer (Default = FALSE)
#' @param pointLayer Adds a point layer (Default = True)
#' @param meanLayer Adds a mean layer (Default = True)
#' @param xlab X axis label (defaults to groupCol)
#' @param ylab Y axis label (defaults to valueCol)
#' @param title Plot title (optional)
#' @param boxColor Color for the boxplot layer (Default = "grey30")
#' @param boxFill Fill Color for the boxplot layer (Default = "deepskyblue3")
#' @param boxAlpha Transparency for the box layer (Default = 0.5)
#' @param violinColor Color for the violin layer (Default = "grey30")
#' @param violinFill Fill Color for the violin (Default = "yellow")
#' @param violinAlpha Transparency for the box layer (Default = 0.5)
#' @param pointColor Color for the point layer (Default = "grey30")
#' @param pointFill Fill color for the point layer (Default = "dodgerblue4")
#' @param pointShape Shape for the point layer (Default = 21; fillable circle)
#' @param pointAlpha Transparency for the box layer (Default = 1)
#' @param boxNotch turn on/off confidence interval notches on boxplots (Default = FALSE)
#' @param boxNotchWidth Set the width of boxnotches (0-1) (Default = 0.8)
#' @param pointSize Size of the points (Default = 4)
#' @param pointJitter Amount to jitter the points (Default = 0) Try 0.2 if you
#'   have alot of points.
#' @param meanColor Color for the mean layer (Default = "red2")
#' @param meanFill Fill color for the mean layer (Default = "goldenrod1")
#' @param meanShape Shape for the mean layer (Default = 21; fillable circle)
#' @param meanAlpha Transparency for the mean layer (Default = 0.7)
#' @param meanSize Size of the mean points (Default = 3)
#' @param baseFontSize The smallest size font in the figure in points. (Default =
#'   12)
#' @param themeStyle "bw" or "grey" which correspond to theme_bw or theme_grey
#'   respectively. Default = bw"
#' @param facet Specifies whether to facet (TRUE) or print individual plots
#'   (FALSE)  (Default = TRUE)
#' @param facetRow Explicitly set the number of rows for the facet plot.
#'   Default behavior will automatically set the columns. (Default = NULL)
#' @param xAngle Angle to set the sample labels on the Xaxis. (Default =  30; Range = 0-90)
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
#'   #simulate some data with row and colnames
#'   groups <- paste("group", factor(rep(1:4, each=100)), sep="")
#'   x <- matrix(rnorm(2400, mean=10), ncol=length(groups))
#'   colnames(x) <- paste("sample", 1:ncol(x), sep="")
#'   rownames(x) <- paste("gene", 1:nrow(x), sep="")
#'   #reformat into tidy dataframe
#'   tidyInt <- tidyIntensity(x,
#'
#'  #Or get data from a DGEobj with RNA-Seq data
#'  log2cpm <- convertCount(dgeObj$counts, unit="cpm", log=TRUE, normalize="tmm")
#'  tidyInt <- tidayIntensity(log2cpm,
#'                            rowidColname="GeneID",
#'                            keyColname="Sample",
#'                            valueColname="Log2CPM")
#'
#'   #Facetted boxplot
#'   obsPlot2(tidyInt, plotByCol="GeneID",
#'            groupCol = "group",
#'            valueCol ="Log2CPM",
#'            pointJitter = 0.1,
#'            facetRow = 2)
#'
#'   #Facetted violin plot
#'   obsPlot2(tidyInt, plotByCol = "GeneID",
#'            violinLayer = TRUE,
#'            boxLayer = FALSE,
#'            groupCo l="group",
#'            valueCol = "Log2CPM",
#'            pointJitter = 0.1,
#'            facetRow = 2)
#'
#'   #Return a list of ggplots for each individual gene
#'   myplots <- obsPlot2(tidyInt, plotByCol="GeneID",
#'            groupCol = "group",
#'            valueCol ="Log2CPM",
#'            pointJitter = 0.1,
#'            facet = FALSE)
#'   #plot one from the list
#'   myplots[[2]]
#'
#' @import ggplot2 magrittr
#' @importFrom dplyr left_join
#' @importFrom assertthat assert_that
#'
#' @export
obsPlot2 <- function(data,
                      plotByCol,
                      groupCol,
                      valueCol,
                      groupOrder=unique(as.character(data[groupCol,, drop=TRUE])),
                      boxLayer = TRUE,
                      violinLayer = FALSE,
                      pointLayer = TRUE,
                      meanLayer = TRUE,
                      xlab=groupCol, ylab=valueCol, title,
                      boxColor = "grey30",
                      boxFill = "deepskyblue3",
                      boxAlpha = 0.5,
                      boxNotch = FALSE,
                      boxNotchWidth = 0.8,
                      violinColor = "grey30",
                      violinFill = "goldenrod1",
                      ViolinAlpha = 0.5,
                      pointColor = "grey30",
                      pointFill = "dodgerblue4",
                      pointShape = 21, #fillable circle
                      pointAlpha = 1,
                      pointSize = 2,
                      pointJitter = 0,
                      meanColor = "red2",
                      meanFill = "goldenrod1",
                      meanShape = 22, #fillable square
                      meanAlpha = 0.7,
                      meanSize = 3,
                      legenPosition = "right",
                      baseFontSize = 12,
                      themeStyle = "grey",
                      facet = TRUE,
                      facetRow,
                      xAngle = 30,
                      scales = "free_y",
                      debug = FALSE
                      )
{

  .addGeoms <- function(MyPlot)
    #note uses global values except for MyPlot
  {

    if (boxLayer==TRUE){
      MyPlot = MyPlot + geom_boxplot(alpha=boxAlpha,
                                     color=boxColor,
                                     fill=boxFill,
                                     notch=boxNotch,
                                     notchwidth = boxNotchWidth,
                                     outlier.shape = outlier.shape,
                                     outlier.size = outlier.size
                                )
    }

    if (violinLayer==TRUE){
      MyPlot = MyPlot + geom_violin(alpha=ViolinAlpha,
                                    color = violinColor,
                                    fill = violinFill)
    }

    if (pointLayer==TRUE){
      if (pointJitter > 0) {
        MyPlot <- MyPlot + geom_point(position = position_jitter(width = pointJitter),
                                    alpha=pointAlpha,
                                    color=pointColor,
                                    fill=pointFill,
                                    size=pointSize,
                                    shape = pointShape)
      } else {
        MyPlot <- MyPlot + geom_point(
                                      alpha=pointAlpha,
                                      color=pointColor,
                                      fill=pointFill,
                                      size=pointSize,
                                      shape=pointShape)
      }
    }

    if (meanLayer==TRUE){
      MyPlot <- MyPlot +
        stat_summary(fun.y=mean, geom="point", shape=meanShape, size=meanSize,
                   color="red", fill = "goldenrod1", alpha=meanAlpha)
    }

    return(MyPlot)
  }

  ### Argument checks
  ###
  assertthat::assert_that(!missing(data),
                          class(data)[[1]] == "data.frame",
                          plotByCol %in% colnames(data),
                          groupCol %in% colnames(data),
                          valueCol %in% colnames(data),
                          all(groupOrder %in% as.character(data[groupCol,,drop=TRUE]))
  )

  #reduce box outliers to a dot if geom_points turned on.
  outlier.size <- 1.5
  outlier.shape <- 19
  if (pointLayer) {
    outlier.size <- 1
    outlier.shape <- "."
  }


  if (debug ==TRUE) browser()



### Plot code here
  if (facet) {

    #set facet columns to sqrt of unique observations (rounded up)
    if (missing(facetRow)) {
      numcol <- data[plotByCol] %>% unique %>% length %>% sqrt %>% ceiling
    } else {
      numcol = facetRow
    }

    # if (numcol > 6) {
    #   warning ("You're putting a lot of plots into a Facet Plot")
    # }

    MyPlot <- ggplot2::ggplot (data, aes_string(x=groupCol, y=valueCol))
    MyPlot <- .addGeoms(MyPlot)
    facetFormula <- stringr::str_c("~", plotByCol, sep=" ")
    MyPlot <- MyPlot + ggplot2::facet_wrap(facetFormula, nrow=numcol, scales=scales)

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

  } else { #individual plots for each Gene returned in a list

      plotlist <- list()

      for (obs in unique(data[[plotByCol]])) {  #for each gene

        dat <- data[data[[plotByCol]] == obs, ] #pull data for one gene
        aplot <- ggplot(dat, aes(x=groupCol, y=valueCol)) + #Samples vs Log2CPM
          xlab(xlab) +
          ylab(ylab) +
          ggtitle(obs) +
          theme_grey() + facetTheme(baseFontSize)
        aplot <- .addGeoms(aplot)
        if (meanLayer == TRUE){
          aplot <- aplot +
            stat_summary(fun.y=mean, geom="point", shape=meanShape, size=meanSize,
                         color="red", fill = "goldenrod1", alpha=1.0)
        }
        #rotate xaxis group labels
        if (xAngle > 0){
          aplot <- aplot + theme(axis.text.x = element_text(angle = xAngle, hjust = 1))
        }
        plotlist[[obs]] <- aplot

      }

      MyPlot = plotlist

  }

  return(MyPlot)
}

