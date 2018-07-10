### Function GOIplot ###
#'
#' Provide a gene list of the genes of interest and return a plot of the LogFC data for each gene.
#' If >1 gene requested, the plot will be faceted which places an upper limit of around 25 or so
#' on the number of genes that can be displayed.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords compare plot ggplot2 logratio scatterplot
#'
#' @param dgeObj A DGEobj workflow datastrucutre containing contrast (topTable) results (required)
#' @param genelist List of ensembl gene IDs for the genes of interest (required)
#' @param plotType One of "point", "bar" or "line" (Default = "bar")
#' @param logRatioCol name of the LogRatio column (Default = "logFC")
#' @param logIntCol name of the LogIntensity column (Default = "AveExpr")
#' @param pvalCol name of the pvalue or FDR column (Default = "P.Value")
#' @param xlab X axis label (default the LogIntensity column name)
#' @param ylab Y axis label (default the LogRatio column name)
#' @param title Plot title (optional)
#' @param symbolSize Size of symbols for Up, no change and Down. default = c(4,
#'   3.99, 4); Note: you cannot choose the exact same size for all three.  But
#'   you can use decimal values.
#' @param symbolShape Shape of the symbols for Up, no change and Down; Default =
#'   c(21, 1, 21) (1 = open circle, 21 = fillable open circle); Note: you cannot
#'   choose the same symbol shape for all three. See
#'   \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Up, NoChange, Down); default = c("black", "grey25",
#'   "grey0") See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#'   Note: You cannot duplicate colors.
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25
#'   are fillable. This will have no effect on other symbols. Default =
#'   c("red3", "grey25", "deepskyblue4") Note: You cannot duplicate colors.
#' @param alpha Controls the transparency of the plotted points (range: 0-1;
#'   default = 0.7)
#' @param zeroRefLine Draw a referenceline at logFC=0 (Default = TRUE)
#' @param refLineThickness Size for the reference line (Default = 0.5)
#' @param baseFontSize The smallest size font in the figure in points. Default =
#'   12
#' @param themeStyle "bw" or "grey" which correspond to theme_bw or theme_grey
#'   respectively. Default = bw"
#' @param footnote optional string placed right justified at bottom of plot.
#' @param footnoteSize applies to footnote. (default = 3)
#' @param footnoteColor applies to footnote. (default = "grey90")
#' @param footnoteJust Value 0-1. 0 is left justified, 1 is right justified, 0.5 is centered. (default=1)
#' @param faceted Set to FALSE to return a list of ggplot objects, one per gene, instead of a single
#'    faceted plot.
#'
#' @return ggplot object or list of ggplot objects
#'
#' @examples
#'    pending
#'
#'
#' @import ggplot2 magrittr dplyr tibble DGEobj
#'
#' @export
GOIplot <- function(dgeObj,
                        logRatioCol = "logFC",
                        logIntCol = "AveExpr",
                        pvalCol = "P.Value",
                        pthreshold=0.01,
                        xlab, ylab, title,
                        symbolSize = 2,
                        symbolShape = 21,
                        symbolColor = "black",
                        symbolFill = "dodgerblue3",
                        alpha = 0.5,
                        zeroRefLine = TRUE,
                        refLineThickness = 1,
                        baseFontSize = 12,
                        themeStyle = "grey",
                        footnote,
                        footnoteSize=3,
                        footnoteColor="black",
                        footnoteJust=1
                        ) {


  #argument checks
  #
  # Make sure specified columns exist
  if (!logRatioCol %in% colnames(df)) {
    stop("LogRatio Ccolumn not found.")
  }
  if (!logIntCol %in% colnames(df)) {
    stop("LogIntensity column not found.")
  }
  if (!pvalCol %in% colnames(df)) {
    stop("Significance measure column not found.")
  }

  #symbol parameters must all be length=3
  if (!length(symbolSize)==3 || !length(symbolShape)==3 ||
      !length(symbolColor)==3 || !length(symbolFill)==3){
    stop("symbol arguments must be length=3 (symbolSize, SymbolShape, SymbolColor, symbolFill)")
  }

#   #all symbolSize values cannot be the same or ggplot throws an error.
#   if (symbolSize[1] == symbolSize[2]) {
#     symbolSize[2] = symbolSize[2] - 0.01
#   }
#
#   #Check for unique colors (all three colors must be unique for the legend to work properly,
#   #even if the fill color is being applied to a non-fillable symbol
#   if (length(unique(symbolColor)) < 3){
#     stop("All 3 symbol colors must be unique")
#   }
#   if (length(unique(symbolFill)) < 3){
#     stop("All 3 symbol fill colors must be unique")
#   }
#   if (length(unique(symbolShape)) < 2){
#     stop("You need to use at least 2 different symbol shapes")
#   }

  if (sizeByIntensity==TRUE){
    #create a column to support sizeByIntensity
    df$LogInt = df[[logIntCol]]
    #set a floor and a ceiling
    df$LogInt[df$LogInt<0] = 0
    df$LogInt[df$LogInt>10] = 10
  }

  names(symbolShape) = c("Increased", "No Change", "Decreased")
  names(symbolSize) = c("Increased", "No Change", "Decreased")
  names(symbolColor) = c("Increased", "No Change", "Decreased")
  names(symbolFill) = c("Increased", "No Change", "Decreased")

  ssc = data.frame(group = c("Increased", "No Change", "Decreased"),
                   symbolShape = symbolShape,
                   symbolSize = symbolSize,
                   symbolColor = symbolColor,
                   symbolFill = symbolFill,
                   order = c(1,3,2),
                   stringsAsFactors = FALSE) %>% arrange(order)
                  #setting the order defines the legend order

  #columns to plot
  #capture the labels from the colname
  xlabel = logRatioCol
  ylabel = paste ("-log10(", pvalCol, ")", sep="")
  #now make the columnames suitable for use with aes_string
  x = make.names(colnames(df)[colnames(df)==logRatioCol])
  colnames(df)[colnames(df)==logRatioCol] = make.names(colnames(df)[colnames(df)==logRatioCol])
  #make a log10significance column and make that the y column
  df$NegativeLogP = -log10(df[,pvalCol])
  y = "NegativeLogP"

  ### DELUXE PLOT: plot groups in different colors/shapes

  #Let's plot the subsets
  DEup = df[[pvalCol]] <= pthreshold & df[[logRatioCol]] > 0
  DEdn = df[[pvalCol]] <= pthreshold & df[[logRatioCol]] < 0
  DEnot = !DEup & !DEdn
  #create group factor column in df
  df$group=NA
  df$group[DEup] = "Increased"
  df$group[DEdn] = "Decreased"
  df$group[DEnot] = "No Change"
  df %<>% left_join(ssc)
  df$group %<>% factor(levels=c("Increased", "Decreased", "No Change"))

  #set an order field to force plotting of NoChange first
  df$order = NA
  df$order[DEup] = 1
  df$order[DEdn] = 1
  df$order[DEnot] = 0

  volcanoPlot <- ggplot (df, aes_string(x=x, y=y)) +
                    aes(shape=group, size=group,
                      color=group, fill=group,
                      order=order) +
      #scale lines tell it to use the actual values, not treat them as factors
      scale_shape_manual(name="Group", guide="legend", labels=ssc$group,
                           values=ssc$symbolShape) +
      scale_size_manual(name="Group", guide="legend", labels=ssc$group,
                          values=ssc$symbolSize) +
      scale_color_manual(name="Group", guide="legend", labels=ssc$group,
                           values=ssc$symbolColor) +
      scale_fill_manual(name="Group", guide="legend", labels=ssc$group,
                          values=ssc$symbolFill) +
      geom_point(alpha=alpha) +
      #box around the legend
      theme(legend.background = element_rect(fill="gray95", size=.5, linetype="dotted"))
###end DELUXE PLOT

  ### Optional Decorations

  if (!is.null(rugColor)){
    volcanoPlot = volcanoPlot + geom_rug(data=df, inherit.aes=FALSE,
                                         color = rugColor,
                                         alpha=rugAlpha,
                                         show.legend=FALSE,
                                         aes_string(x=x, y=y))
  }

  if (sizeByIntensity==TRUE) {
    volcanoPlot = volcanoPlot + aes(size=LogInt) +
      scale_size_continuous()
  }

  if (!is.null(pthresholdLine)){
    volcanoPlot = volcanoPlot +
      geom_hline(yintercept=-log10(pthreshold), color=pthresholdLine,
                 alpha=0.5, size=refLineThickness)
  }

  if (!is.null(foldChangeLines)){
    volcanoPlot = volcanoPlot +
      geom_vline(xintercept=foldChangeLines, color=symbolFill["Increased"],
                 alpha=0.5, size=refLineThickness) +
      geom_vline(xintercept=-foldChangeLines, color=symbolFill["Decreased"],
                 alpha=0.5, size=refLineThickness)
  }

  ### Add Labels

  if (is.null(xlab)){ #use colname unless supplied as argument
    volcanoPlot = volcanoPlot + xlab(xlabel)
  } else {
    volcanoPlot = volcanoPlot + xlab(xlab)
  }
  if (is.null(ylab)){
    volcanoPlot = volcanoPlot + ylab(ylabel)
  } else {
    volcanoPlot = volcanoPlot + ylab(ylab)
  }
  if (!is.null(title)){
    volcanoPlot = volcanoPlot +
      ggtitle(title)
  }

  #Set the font size before placing the legend
  if (tolower(themeStyle) == "bw") {
    volcanoPlot = volcanoPlot + theme_bw() + baseTheme(baseFontSize)
  } else {
    volcanoPlot = volcanoPlot + theme_grey() + baseTheme(baseFontSize)
  }

  volcanoPlot = setLegendPosition(volcanoPlot, legendPosition, themeStyle)

  #footnote
  if (!missing(footnote))
    volcanoPlot <- volcanoPlot (CompPlot, footnoteText=footnote,
                                footnoteSize=footnoteSize,
                                footnoteColor="black",
                                footnoteJust=footnoteJust)

  return(volcanoPlot)
}

