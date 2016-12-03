### Function ggplotMDS ###
#' Function ggplotMDS
#'
#' This is a wrapper around the plotMDS function that generates the plot with
#' ggplot2 instead of base graphics.
#'
#' The underlying plotMDS function uses a default of top=500 to use the top 500
#' highest fold change genes for the analysis.  The method runs quick and I've
#' found situations where a larger value produced a more stable result.  Thus
#' I've change the default to top=Inf.  You should play with this parameter and
#' at least try a number that is close to the number of differential genes
#' detected in your data.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords compare plot ggplot2 logratio scatterplot
#'
#' @param DGElist A DGEList object taken after running edgeR TMM normalization
#' @param top Number of most variant genes to include (default = Inf)
#' @param labels Text labels for the samples.  These should be short
#'   abbreviations of the sample identifiers so that the plot is not too
#'   cluttered.  If NULL, the column names from DGElist are used.
#' @param labelSize Size for the text labels in the plot (default = 2)
#' @param textColor Color for the text labels in the plot (default = "blue2")
#' @param vlineIntercept X intercept of vertical line (Default = Null)
#' @param hlineIntercept Y intercept of horizontal line (Default = Null)
#' @param reflineColor Color for the horizontal and vertiacl reference lines
#'   (Default = "darkgoldenrod1")
#' @param reflineSize Thickness of the reference lines (Default = 0.5)
#' @param baseFontSize Base fontsize for the plot (Default = 12)
#'
#' @param themeStyle One of "grey" or "bw" (Default = "grey")
#' @param ... arguments from plotMDS.DGEList (see ?plotMDS.DGEList in package
#'    edgeR)
#'
#' @return A list with two elements, the ggplot object and the mds object returned
#'    by the plotMDS function.
#'
#' @examples
#'
#'      #Plot the first two dimensions using all genes
#'      MyMDS = ggplotMDS (MyDGEList)
#'      #display the plot and save a png file.
#'      printAndSave (MyMDS[[1]], file = "MyMDS.PNG", width=5, height=5)
#'
#'      #plot the 2nd and 3rd dimensions using the top 1000 genes
#'      MyMDS = ggplotMDS (MyDGEList, dim.plot=c(2,3) ndim =3)
#'
#' @import ggplot2 magrittr limma edgeR
#'
#' @export
ggplotMDS <- function(DGElist,
                      top = Inf,
                      labels = NULL,
                      labelSize = 2,
                      title = NULL,
                      textColor = "blue2",
                      hlineIntercept = NULL,
                      vlineIntercept = NULL,
                      reflineColor = "darkgoldenrod1",
                      reflineSize = 0.5,
                      baseFontSize = 12,
                      themeStyle = "grey",
                      ...
                        )
{

  #argument checks
  if (!tolower(class(DGElist)[[1]]) == "dgelist") {
    stop("DGElist must be class DGEList. Try yourSLOA$DGElist.")
  }

  #defaults
  if (!exists("dim.plot")){
    dim.plot <- c(1,2)
  }
  if (!exists("ndim")){
    ndim <- 2
  }
  if (!exists("gene.selection")){
    gene.selection <- "pairwise"
  }
  if (!exists("Xlab")){
    Xlab <- NULL
  }
  if (!exists("Ylab")){
    Ylab <- NULL
  }
  if (!exists("pch")){
    pch <- NULL
  }
  if (!exists("cex")){
    cex <- 1
  }
  if(!exists("method")){
    method <- "logFC"
  }
  if (!exists("prior.count")){
    prior.count <- 2
  }
  if (is.null(labels)){
    labels <- colnames(DGElist)
  }
  if (is.null(title)){
    title <- "MDS Plot"
  }

  pdf(NULL) #suppress the plot and just capture the output
  mds <- plotMDS(DGElist, top = top, labels = labels, pch = pch,
                 cex = cex, dim.plot = dim.plot, ndim = ndim,
                 gene.selection = gene.selection,
                 xlab = Xlab, ylab = Ylab)
  invisible(dev.off())

  #draw the plot
  xydat = cbind(x=mds$x, y=mds$y) %>% as.data.frame
  if (is.null(Xlab) && is.null(Ylab)){
    xylab = paste(mds$axislabel, mds$dim.plot, sep=" ")
  } else {
    xylab = c(Xlab, Ylab)
  }

  mdsplot = ggplot(xydat, aes(x=x, y=y)) +
    #geom_point(color="darkblue", shape=21, fill = "blue") +
    geom_text(color=textColor, label=labels, size=labelSize) +
    xlab (xylab[1]) +
    ylab (xylab[2]) +
    ggtitle (title)

  #place an annotation on the bottom left of the plot
  xrange = getXrange(mdsplot)
  yrange = getYrange(mdsplot)
  #put the annotation 10% from xmin
  xpos = xrange[1] + ((xrange[2] - xrange[1]) * 0.1 )
  alabel = paste("top ", mds$top, " genes : gene.selection = ",
                 mds$gene.selection, sep="")
  mdsplot = mdsplot + annotate ("text", x = xpos, y = yrange[1],
                                label = alabel, hjust=0,
                                size=rel(2.5), color="grey30")

  if (!is.null(hlineIntercept)){
    mdsplot <- mdsplot + geom_hline (yintercept = hlineIntercept,
                                     color = reflineColor,
                                     size=reflineSize)
  }
  if (!is.null(vlineIntercept)){
    mdsplot <- mdsplot + geom_vline (xintercept = vlineIntercept,
                                     color = reflineColor,
                                     size=reflineSize)
  }

  if (tolower(themeStyle) %in% c("grey", "gray")){
    mdsplot <- mdsplot + greyTheme(baseFontSize)
  }
  #mdsplot
  MDS <- list(plot=mdsplot, mdsobj=mds)

  return(MDS)
}

