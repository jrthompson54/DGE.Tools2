### Function cdfPlot ###
#' Deluxe CDF Plots
#'
#' CDF plots are a good complement to pvalue histograms as a way to evaluate
#' model performance and examine support for differential expression. Results
#' are ranked by pvalue on the x-axis and the pvalue plotted on the y-axis.
#' Since pvalue distributions should be flat, this type of plot should produce a
#' straight line.  Any observations that fail to meet the null hypothesis will
#' appear as a break in the line at the low end of the curve.
#'
#' This function is designed to take topTable dataframes and display the
#' corresponding CDF plot. Data for the pvalues below 0.1 (user settable via
#' pvalMax argument) are show in a full size plot. An inset figure shows the
#' whole pvalue scale and highlights the portion shown in the full plot.  Points
#' below 0.01 are a different color by default (threshold set by pthreshold
#' argument; shape/color attributes customizable through arguments).
#'
#' Two arguments control the output.  printPlot = TRUE outputs the compound plot
#' to the console/knitr.  Provide a filename with a .PNG extension to save the
#' compound plot to a PNG file.  Font sizes are increase 1.5 fold for the PNG
#' file. The default baseFontSize = 12 produces suitable output for knitr PDFs
#' and the 1.5 fold adjustment produces a PNG file suitable for PPT
#' presentations.
#'
#' \strong{Data Structure for the input dataframe:}
#'
#' The defaults are set for dataframes produced by topTable.  The column
#' "P.Value" is used by default to accommodate the column names used in topTable
#' dataframes.  Any other dataframe can be used with by explicitly defining the
#' pvalue column name with the appropriate argument.
#'
#' Sensible defaults are chosen for symbols (Size, Shape, Color and Fill).
#' Optionally, there are arguments to allow you to fiddle with these settings.
#' A length of 2 is required for these arguments which applies the attributes in
#' this order: Significant, Not Significant.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords compare plot ggplot2 logratio scatterplot
#'
#' @param df A dataframe with LogRatio and LogIntensity columns and optionally a pvalue or FDR column.
#' @param pvalCol name of the pvalue or FDR column (Default = "P.Value")
#' @param pvalMax Limit the range of the main plot (Default = 0.10)
#' @param pthreshold Used to color points (default = 0.01)
#' @param xlab X axis label (default is "LogIntensity column name"Rank")
#' @param ylab Y axis label (default is Pvalue column name)
#' @param title Plot title (optional)
#' @param insetTitle Title for the inset plot (Optional)
#' @param symbolSize Size of symbols for Up, no change and Down. default = c(4, 1, 4);
#'        Note: you cannot choose the exact same size for all three.  But you can use decimal
#'        values.
#' @param symbolShape Shape of the symbols for Up, no change and Down; Default = c(21, 20, 21)
#'        (20 = filled circle, 21 = fillable open circle); Note: you cannot choose the same
#'        symbol shape for all three.
#'        See \url{http://www.cookbook-r.com/Graphs/Shapes_and_line_types}
#' @param symbolColor c(Up, NoChange, Down); default = c("black", "grey25", "grey0")
#'        See \url{http://research.stowers-institute.org/efg/R/Color/Chart}
#'        Note: You cannot duplicate colors.
#' @param symbolFill Set the fill color for the symbols. Note only symbols 21-25 are fillable.
#'        This will have no effect on other symbols.
#'        Default = c("red3", "grey25", "deepskyblue4")
#'        Note: You cannot duplicate colors.
#' @param alpha Controls the transparency of the plotted points (0-1; default =
#'   0.5)
#' @param referenceLine Color for an horizontal line drawn at the pthreshold
#'   (Default = NULL; NULL disables, set to desired color to enable)
#' @param refLineThickness Set thicknesss of the reference line (Default=1)
#' @param baseFontSize The smallest size font in the figure in points. (Default
#'   = 12)
#' @param themeStyle "bw" or "grey" which correspond to theme_bw or theme_grey
#'   respectively. Default = "grey"
#' @param printPlot Specify printing the combined plot to the console/knitr
#'   (Default = TRUE)
#' @param plotFile Provide a filename with .PNG extension to save the an image
#'   file. Font size doubled for PNG file output.
#' @param footnote optional string placed right justified at bottom of plot.
#' @param footnoteSize applies to footnote. (default = 3)
#' @param footnoteColor applies to footnote. (default = "black")
#' @param footnoteJust Value 0-1. 0 is left justified, 1 is right justified, 0.5 is centered. (default=1)
#'
#' @return A list containing mainplot, insetplot and viewport. You can then
#'   reconstruct the plot with: print(MyList$main); print(inset,
#'   vp=Mylist$viewport)
#'
#' @examples
#'    Plot to console and PNG file (df is a topTable dataframe)
#'
#'    cdfPlot(df, title = "My CDF Plot", printFile = "MyCDFplot.PNG")
#'
#' @import ggplot2 magrittr dplyr grid
#'
#' @export
cdfPlot <- function(df,
                        pvalCol = "P.Value",
                        pthreshold=0.01,
                        xlab=NULL, ylab=NULL, title=NULL, insetTitle=NULL,
                        symbolSize = c(2, 1),
                        symbolShape = c(20, 20),
                        symbolColor = c("red3", "deepskyblue4"),
                        symbolFill = c("red3", "deepskyblue4"),
                        alpha = 1,
                        referenceLine = NULL,
                        refLineThickness = 1,
                        legendPosition = "se",
                        baseFontSize = 12,
                        viewportX = 0.15,
                        viewportY = 0.93,
                        viewportWidth = 0.35,
                        themeStyle = "grey",
                        pvalMax = 0.10,
                        printPlot = TRUE,
                        plotFile,
                        footnote,
                        footnoteSize=3,
                        footnoteColor="black",
                        footnoteJust=1
                        )
{

  #argument checks
  #
  # Make sure specified columns exist

  if (!pvalCol %in% colnames(df)) {
    stop("Significance measure column not found.")
  }

  #symbol parameters must all be length=2
  if (!length(symbolSize)==2 || !length(symbolShape)==2 ||
      !length(symbolColor)==2 || !length(symbolFill)==2){
    stop("symbol arguments must be length=2 (symbolSize, SymbolShape, SymbolColor, symbolFill)")
  }

  names(symbolShape) <- c("Significant", "Not Significant")
  names(symbolSize) <- c("Significant", "Not Significant")
  names(symbolColor) <- c("Significant", "Not Significant")
  names(symbolFill) <- c("Significant", "Not Significant")

  ssc <- data.frame(group = c("Significant", "Not Significant"),
                   symbolShape = symbolShape,
                   symbolSize = symbolSize,
                   symbolColor = symbolColor,
                   symbolFill = symbolFill,
                   order = c(1,2),
                   stringsAsFactors = FALSE) %>% arrange(order)
                    #setting the order defines the legend order

  #columns to plot
  #capture the labels from the colname
  xlabel <- "Rank"
  ylabel <- pvalCol
  x <- xlabel
  y <- ylabel

  ### Combo PLOT: full data inset, most significant data in main plot

  #rank by pvalue
  df %<>% dplyr::arrange_(pvalCol)
  df$Rank <- c(1:nrow(df))

  #Let's plot the pvalue subsets
  Sig <- df[[pvalCol]] <= pthreshold
  NotSig <- df[[pvalCol]] > pthreshold

  #create group factor column in df  #needed for legend I think
  df$group <- NA
  df$group[Sig] <- "Significant"
  df$group[NotSig] <- "Not Significant"
  df %<>% left_join(ssc)
  df$group %<>% factor(levels=c("Significant", "Not Significant"))

  #set an order field to force plotting of NotSig first
  df$order <- NA
  df$order[NotSig] <- 1
  df$order[Sig] <- 2

  #rows to include in the zoomed in plot
  #subsetRows = nrow(df) %>% multiply_by(subset) %>% round
  subsetRows <- sum(df[[y]] <= pvalMax)
  dfsubset <- df[1:subsetRows,]

  #plot subset percent of the data for the main plot
  cdfMain <- ggplot (dfsubset, aes_string(x=x, y=y)) +
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
    geom_point(alpha=alpha)

  ### Optional Decorations
  if (!is.null(referenceLine)){
    cdfMain <- cdfMain +
      geom_hline(yintercept=pthreshold, color=referenceLine,
                 size=refLineThickness, alpha=0.5)
  }

  ### Add Labels

  if (is.null(xlab)){ #use colname unless supplied as argument
    cdfMain <- cdfMain + xlab(xlabel)
  } else {
    cdfMain <- cdfMain + xlab(xlab)
  }
  if (is.null(ylab)){
    cdfMain <- cdfMain + ylab(ylabel)
  } else {
    cdfMain <- cdfMain + ylab(ylab)
  }
  if (!is.null(title)){
    cdfMain <- cdfMain +
      ggtitle(title)
  }

  #Set the font size before placing the legend
  if (tolower(themeStyle) == "bw") {
    cdfMain <- cdfMain + theme_bw(baseFontSize) #+ baseTheme(baseFontSize)
  } else {
    cdfMain <- cdfMain + theme_grey(baseFontSize) #+ baseTheme(baseFontSize)
  }

  cdfMain <- setLegendPosition(cdfMain, legendPosition, themeStyle)

  #footnote
  if (!missing(footnote))
    cdfMain <- addFootnote (cdfMain, footnoteText=footnote,
                            footnoteSize=footnoteSize,
                            footnoteColor="black",
                            footnoteJust=footnoteJust)

  ###end Main PLOT

  ### Set up the inset plot with All Data
  cdfInset <- ggplot (df, aes_string(x=x, y=y)) +
    aes(shape=group, size=group,
        color=group, fill=group,
        order=order) +
    #scale lines tell it to use the actual values, not treat them as factors
    scale_shape_manual(name="Group", guide="none", labels=ssc$group,
                       values=ssc$symbolShape) +
    scale_size_manual(name="Group", guide="none", labels=ssc$group,
                      values=ssc$symbolSize) +
    scale_color_manual(name="Group", guide="none", labels=ssc$group,
                       values=ssc$symbolColor) +
    scale_fill_manual(name="Group", guide="none", labels=ssc$group,
                      values=ssc$symbolFill) +
    geom_rect(xmin = 0, xmax = subsetRows,
              ymin = 0, ymax = max(dfsubset[[y]]), color = "lightblue",
              fill = "lightblue", alpha=0.2) +
    geom_point(alpha=alpha)

  ### Add Labels

  if (is.null(xlab)){ #use colname unless supplied as argument
    cdfInset <- cdfInset + xlab(xlabel)
  } else {
    cdfInset <- cdfInset + xlab(xlab)
  }
  if (is.null(ylab)){
    cdfInset <- cdfInset + ylab(ylabel)
  } else {
    cdfInset <- cdfInset + ylab(ylab)
  }
  if (!is.null(insetTitle)){
    cdfInset <- cdfInset +
      ggtitle(insetTitle)
  }

  #Adjust font size for inset.
  factor <- 4/baseFontSize


  #Set the font size
  if (tolower(themeStyle) == "bw") {
    cdfInset <- cdfInset + theme_bw() + baseTheme(baseFontSize*factor)
  } else {
    cdfInset <- cdfInset + theme_grey() +baseTheme(baseFontSize*factor)
  }

  #Now arrange plot inside a plot

  #Adjust viewport Y if main Title present
  vy <- viewportY
  if (!is.null(title)){
    adjust <- (baseFontSize/150)
    vy <- viewportY - adjust
    #viewportX %<>% add(adjust)
  }
  #A viewport taking up a fraction of the plot area (upper left)
  vp <- viewport(width = viewportWidth, height = viewportWidth,
                 x = viewportX, y = vy,
                 just = c("left", "top"))

  if (printPlot == TRUE) {
    print(cdfMain)
    print(cdfInset, vp=vp)
  }

  ### PNG file output
  if (!missing(plotFile)){

    #increase font size for PNG file
    pngFontSize <- baseFontSize * 1.5

    #Adjust viewport X for larger font
    adjust <- (pngFontSize/480)
    viewportX %<>% add(adjust)

    #Adjust viewport Y if main Title present
    vy <- viewportY
    if (!is.null(title)){
      adjust <- (baseFontSize/100)
      vy <- viewportY - adjust
    }

    #Shrink the legend by 50%
    cdfMain <- cdfMain +
      theme(legend.key.size = unit((7.5/pngFontSize), "lines"),
            legend.text = element_text(size = rel(6/pngFontSize)),
            legend.title = element_text(size = rel((7/pngFontSize)), face = "bold")
      )

    #A viewport taking up a fraction of the plot area (upper left)
    vp <- viewport(width = viewportWidth, height = viewportWidth,
                   x = viewportX, y = vy,
                   just = c("left", "top"))
    png(file=plotFile, width=6, height=4, units = "in", res = 300)
    print(cdfMain + baseFont(pngFontSize))
    print(cdfInset + baseFont(pngFontSize/3), vp=vp)
    invisible ( dev.off() )
  }

#   ### create a function to print the combined plot
#   comboPlot <- function(baseFontSize=12)
#   {
#     #Adjust viewport ...
#     adjust = (baseFontSize/480)
#     if (!is.null(title)){
#       viewportY %<>% subtract(adjust)
#     }
#     #adjust X for fontsize
#     viewportX %<>% add(adjust)
#
#     #A viewport taking up a fraction of the plot area (upper left)
#     vp <- viewport(width = viewportWidth, height = viewportWidth,
#                    x = viewportX, y = viewportY,
#                    just = c("left", "top"))
#     print(cdfMain + baseFont(baseFontSize))
#     smallFontSize <- baseFontSize %>% multiply_by(0.5)
#     print(cdfInset + baseFont(smallFontSize), vp=vp)
#   }

  MyList = c(main=cdfMain, inset=cdfInset, viewport=vp)

  return(MyList)
}

