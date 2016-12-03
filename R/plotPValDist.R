### Function plotPvalDist ###
#' Function  plotPvalDist
#'
#' Generate a facet plot (or optionally individual plots) from a dataframe of
#' numbers.  Intended to perform histogram analysis of pvalue distributions.
#' But should be useful for any dataframe of numeric columns.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID
#'
#' @param P.Val a matrix or dataframe of numeric data; col=samples
#' @param Facet TRUE adds the samplename as a title bar on each facet
#' @param savePlot  TRUE saves the plot(s) to .PNG files
#' @param fileNames  Only relevant when Facet = FALSE. Default behavior takes png file
#'		names from colnames.  Supply the FileNames parameter with a list of filenames
#'		to override the default naming.  Must be same length as ncol(P.Val)
#'
#' @return A ggplot2 object or a list of plots
#'
#' @examples
#' MyDataframe = plotPvalDist (MyPvalDF, Facet = TRUE, SavePlot = TRUE)
#'
#' @import ggplot2 reshape2
#' @export
plotPvalDist <- function (P.Val, Facet = TRUE, savePlot = TRUE, fileNames = NULL) {

  #JRT 5Aug2015
  #Input:
  #   minimally a pvalue matrix or dataframe with sample names as colnames.
  #   and GeneIDs as rownames.
  #
  #   Facet=FALSE will force individual plots instead of a facet plot
  #
  #   Plots will be saved to png files also unless SavePlot = FALSE.
  #   The facet plot will be names "PValue_Dist_facet.png"
  #
  #   Individual plot filenames will be samplename.png where samplename is taken from
  #   colnames unless you provide a vector of base filenames (without the .png)
  #
  #Output: plot pvalue distributions as a faceted plot or optionally as individual plots
warning ("Warning: plotPvalDist is deprecated and no longer supported.  Please use plotPvalHist instead.")

  mytheme_facet = ggplot2::theme(
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    axis.title.x = element_text(face="bold", colour="Black", size=20),
    axis.title.y = element_text(face="bold", colour="Black", size=20),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    plot.title = element_text(lineheight=.8, face="bold", size = 24),
    legend.text = element_text(colour="Black", size=16),
    legend.title = element_text(colour="Black", size=14),
    strip.text.x = element_text(size=16),
    strip.text.y = element_text(size=16)
  )

  mytheme = ggplot2::theme(axis.text.x = element_text(size=18),
                           axis.text.y = element_text(size=18),
                           axis.title.x = element_text(face="bold", colour="Black", size=20),
                           axis.title.y = element_text(face="bold", colour="Black", size=20),
                           plot.title = element_text(lineheight=.8, face="bold", size = 24),
                           legend.text = element_text(colour="Black", size=16),
                           legend.title = element_text(colour="Black", size=14),
                           strip.text.x = element_text(size=16),
                           strip.text.y = element_text(size=16)
  )


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
    if (is.null(FileNames)) {
      FileNames = "PValue_Dist_facet.png"
    }

    Dist_Pval_Facet <- ggplot2::ggplot(data=P.Val, aes(x=Pval)) +
      ggplot2::geom_density(alpha = 0.6) +
      ggplot2::xlab("Pvalue") +
      ggplot2::ylab("Density") +
      ggtitle("Pvalue Distributions") +
      ggplot2::scale_fill_brewer(palette="Set1") +
      # facet_grid(Levels ~ ., scales="free") +
      ggplot2::facet_wrap(~ Levels, nrow = numrow, scales="free") +
      mytheme_facet


    if (SavePlot) {
      #print(Dist_Pval_Facet)

      png(file=FileNames[1], width=8, height=6, units = 'in', res = 300)
      print(Dist_Pval_Facet)
      invisible ( dev.off() )
    }
    return(Dist_Pval_Facet)

  } else { #Print individual plots
    #decide on filenames
    if (!is.null(FileNames)) {
      #confirm we have the right number of filenames
      if (!length(FileNames) == ncol(P.Val)) {
        stop("Number of FileNames does not match column count!")
      } else {
        FileNames = SampNames
        #substitute underscores for colons
        FileNames = gsub(":", "_", FileNames)
      }
    }

    #run off the plots
    plotlist = list()
    for (i in 1:NumSamples) {
      if (is.null(FileNames))
      f = FileNames[i]
      s = SampNames[i]
      MyPVal = filter (P.Val, grepl(s, Levels))

      Dist_Pval <- ggplot2::ggplot(data=MyPVal, aes(x=Pval)) +
        ggplot2::geom_density(alpha = 0.6) +
        ggplot2::xlab("Pvalue") +
        ggplot2::ylab("Density") +
        ggplot2::ggtitle(paste("Pvalue Distributions\n", s)) +
        #scale_fill_brewer(palette="Set1") +
        mytheme

      print(Dist_Pval)
      plotlist[i] = Dist+Pval

      if (saveplot) {
        png(file=paste(f, ".png", sep=""),width=8,height=6, units = 'in', res = 300)
        print(Dist_Pval)
        invisible ( dev.off() )
      }
      return(plotlist)
    }
  }
}
