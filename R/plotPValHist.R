### Function plotPvalHist ###
#' Function  plotPvalHist
#'
#' Generate a facet plot (or optionally individual plots) from a dataframe of
#' numbers.  Intended to perform histogram analysis of pvalue distributions.
#' But should be useful for any dataframe of numeric columns.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID
#'
#' @param P.Val a matrix or dataframe of numeric data; col=samples
#' @param Facet Set to FALSE to print individual plots instead of a facetted plot. (Default = TRUE)
#' @param savePlot TRUE saves the plot(s) to .PNG files
#' @param fileNames Pass a single filename for a facetted plot to save the
#'        plot to a file.  Pass a list of filenames if Facet=FALSE (length must equal
#'        sample count).
#' @param binWidth Range is always 0-1 for pvalues. So this defaults to 0.02
#' @param baseFontSize Set the base font size for the plot using theme_grey (Default=12)
#' @param facetFontSize control font size on the individual plot headers
#'   (default=NULL). Allows for fine tuning if the theme_grey default doesn't
#'   work for you. 10 seems good for knitr output (dependant on length of your
#'   sample titles); bigger (e.g. 14) works better for PPT.
#' @param alpha Set the transparency. (Default = 0.6)
#' @param color Color for the histogram outline (default = "dodgerblue3")
#' @param fill Fill color for the histogram (default = "dodgerblue3")
#'
#' @return A ggplot2 object if facet=TRUE or a list of plots if facet=FALSE. (Default = TRUE)
#'
#' @examples
#' #Print to console, all defaults
#' plotPvalHist (MyPvalMatrix)
#'
#' #Print just to a PNG with some arguments
#' myplot = plotPvalHist (MyPValMatrix, savePlot = T, fileNames = "MyPlot.PNG",
#'                facetFontSize= 14)
#'
#' #Print the previous plot to the console
#' print(myplot)
#'
#' @import ggplot2 reshape2
#' @export
plotPvalHist <- function (P.Val, Facet = TRUE,
                          savePlot = FALSE,
                          fileNames = NULL,
                          binWidth = 0.02,
                          baseFontSize=12,
                          facetFontSize=NULL,
                          alpha = 0.6,
                          color = "dodgerblue3",
                          fill = "dodgerblue3") {

  #JRT 5Aug2015
  #Input:
  #   minimally a pvalue matrix or dataframe with sample names as colnames.
  #   and GeneIDs as rownames.
  #
  #   Facet=FALSE will force individual plots instead of a facet plot
  #
  #   Plots will be saved to png files also unless savePlot = FALSE.
  #   The facet plot will be names "PValue_Dist_facet.png"
  #
  #   Individual plot filenames will be samplename.png where samplename is taken from
  #   colnames unless you provide a vector of base filenames (without the .png)
  #
  #Output: plot pvalue distributions as a faceted plot or optionally as individual plots

  #9Dec2015  Update to used scaled themes


  if (is.matrix(P.Val)) {
    P.Val %<>% as.data.frame
  }

  NumSamples = ncol(P.Val)
  SampNames = colnames(P.Val)

  #set up Tall format
  P.Val$GeneID = rownames(P.Val)
  P.Val %<>% reshape2::melt(variable.name = "Levels", value.name = "Pval") #(P.Val, id.vars="GeneID")

  if (Facet) {

    numcol = 3
    numrow = NumSamples/numcol %>% ceiling
    if (is.null(fileNames)) {
      fileNames = "PValueHistFacet.png"
    }

    Hist_Pval_Facet <- ggplot2::ggplot(data=P.Val, aes(x=Pval)) +
      ggplot2::geom_histogram(alpha = alpha, fill=fill, color = color,
                              binwidth=binWidth) +
      ggplot2::xlab("Pvalue") +
      ggplot2::ylab("Count") +
      ggtitle("Pvalue Histograms") +
      ggplot2::scale_fill_brewer(palette="Set1") +
      ggplot2::facet_wrap(~ Levels, nrow = numrow, scales="free") +
      theme_grey() + baseTheme(baseFontSize)

      if (!is.null(facetFontSize)) {
       Hist_Pval_Facet <- Hist_Pval_Facet +
        theme(strip.text.x = element_text(size = facetFontSize,
                                        colour = "red", angle = 0))
      }

    if (savePlot) {
      #print(Hist_Pval_Facet)
      png(file=fileNames[1], width=8, height=6, units = 'in', res = 300)
      print(Hist_Pval_Facet)
      invisible ( dev.off() )
    }
    return(Hist_Pval_Facet)

  } else { #Print individual plots

    #decide on filenames
    if (!is.null(fileNames)) {
      #confirm we have the right number of filenames
      if (!length(fileNames) == ncol(P.Val)) {
        stop("Number of fileNames does not match column count!")
      } else {
        fileNames = SampNames
        #substitute underscores for colons
        fileNames = gsub(":", "_", fileNames)
      }
    }

    #run off the plots
    plotlist = list()
    for (i in 1:NumSamples) {
      if (is.null(fileNames))
      f = fileNames[i]
      s = SampNames[i]
      MyPVal = filter (P.Val, grepl(s, Levels))

      Hist_Pval <- ggplot2::ggplot(data=MyPVal, aes(x=Pval)) +
        ggplot2::geom_histogram(alpha = alpha, fill = fill, color = color,
                                binwidth=binWidth) +
        ggplot2::xlab("Pvalue") +
        ggplot2::ylab("Count") +
        ggplot2::ggtitle(paste("Pvalue Histogram\n", s)) +
        theme_grey() + baseTheme(baseFontSize)

      #print(Hist_Pval)
      plotlist[[i]] = Hist_Pval

      if (savePlot) {
        png(file=paste(f, ".png", sep=""),width=8,height=6, units = 'in', res = 300)
        print(Hist_Pval)
        invisible ( dev.off() )
      }
      return(plotlist)
    }
  }
}
