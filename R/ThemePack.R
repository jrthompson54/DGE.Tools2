#JRT 1Dec2015
#
#Theme Pack Overview
#
#you can use baseFont to change just the base font size without disturbing other
#features of you plot in order to get a font size appropriate for the medium.  I
#suggest 12 point for knitr reports and 24 point for powerpoint presentations.
#
#baseTheme() is a basic theme that sets relative fonts sizes as follows:
#  axis scale numbering = 1.0
#  axis titles = 1.25
#  plot title = 1.5
#  legend text = 1.0
#  legend title = 1.2

#The function printAndSave, takes the name of a
#plot object and a filename and prints it to the console with small font
#and saves it to a file with a larger font. For printAndSave to work best,
#you should baseTheme to establish relative font sizes.
#the default font sizes are 12 for printing reports and 24 for saved graphics
#intended for powerpoint.
#
#I also defined a few trivial functions for tweaking my graphics because
#I can never remember the syntax and spend too much time googling to
#look up the syntax on these formatting tasks.  These include:
#
#   NoXlab, NoYlab: turn of axis labels
#   NoLegend: turn of a legend
#   NoXGrid, NoYGrid:  turn off gridlines






### Function baseTheme ###
#' Function  baseTheme
#'
#' A basic theme for individual plots that only sets relative font sizes for
#' common graphic elements.  Apply baseTheme as your last layer.  In particular,
#' theme_bw, theme_grey, will undo some of the changes set by baseTheme.  After
#' applying baseTheme, you can use baseFont to adjust font sizes for different 
#' output purposes (e.g. larger for PPT, smaller for knitr).
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @param base_size Size for the basefont in points (Default = 18)
#' @param base_family Set the font family
#' @param scale_legend Scales Legend.Text to 10/base_size lines (Default=TRUE)
#'
#' @return A ggplot theme that can be added to a plot object
#'
#' @examples
#' Myggplot = myggplot + baseTheme(12)  #start with a 12pt font theme
#' Myggplot = myggplot + baseFont(24)  #change plot to 24pt base.
#'
#' @import ggplot2
#'
#' @export
baseTheme <- function(base_size=18, base_family="", scale_legend = TRUE)
{
  #This is my attempt to automatically scale legend sizes.  I want 1 line at 12 point
  #and 0.5 lines at 24 point.  So let's try 12/base_size as a scaling factor.
  if (scale_legend == TRUE && base_size>17) {
    Legend.ScaledSize <- 14/base_size
    Legend.ScaledFont <- 12/base_size
  } else {
    Legend.ScaledSize = 1.2 #theme_grey defaults
    Legend.ScaledFont = 0.8
  }

	baseTheme <- theme(  #base size is the axis tick mark labels; other elements are scaled
	    axis.text.x = element_text(size=rel(1.0)),
	    axis.text.y = element_text(size=rel(1.0)),
	    axis.title.x = element_text(size=rel(1.25), vjust = 0.5, hjust=0.5, color="black"),
	    axis.title.y = element_text(size=rel(1.25), vjust = 0.5, hjust=0.5, color="black"),
	    plot.title = element_text(face="bold", size = rel(1.5)),
	    legend.text = element_text(colour="Black", size=rel(Legend.ScaledFont)),
	    legend.title = element_text(colour="Black", size=rel(1.2)),
	    legend.title.align = 0.5,
	    legend.key.size = unit(Legend.ScaledSize, "lines"),
	    strip.text.x = element_text(size=rel(1.6)),
	    strip.text.y = element_text(size=rel(1.6)),
	    text = element_text(size=base_size, family = base_family)
	)
}

### Function facetTheme ###
#' Function  facetTheme
#'
#' A basic theme for facetted plots that sets relative font sizes
#' for common graphic elements.  Use in conjunction with baseFont
#' to easily adjust all fonts in a plot.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @param base_size Size for the basefont in points (Default = 18)
#' @param base_family Set the font family
#'
#' @return A ggplot theme that can be added to a plot object
#'
#' @examples
#' Myggplot = myggplot + facetTheme(12)  #start with a 12pt font theme
#' Myggplot = myggplot + baseFont(24)  #change plot to 24pt base.
#'
#' @import ggplot2
#'
#' @export
facetTheme <- function(base_size=18, base_family=""){
	facetTheme = theme(
	    axis.text.x = element_text(size=rel(1.0)),
	    axis.text.y = element_text(size=rel(1.0)),
	    axis.title.x = element_text(face="bold", colour="Black", size=rel(2.0)),
	    axis.title.y = element_text(face="bold", colour="Black", size=rel(2.0)),
	    panel.grid.minor.x=element_blank(),
	    panel.grid.major.x=element_blank(),
	    plot.title = element_text(lineheight=.8, face="bold", size = rel(2.4)),
	    legend.text = element_text(colour="Black", size=rel(1.6)),
	    legend.title = element_text(colour="Black", size=rel(1.4)),
	    legend.title.align = 0.5,
	    strip.text.x = element_text(size=rel(1.6)),
	    strip.text.y = element_text(size=rel(1.6)),
	    text = element_text(size=base_size, family = base_family)
	  )
}

### Function baseFont ###
#' Function  baseFont
#'
#' Sets just the font size and optionally font family without changing other
#' graphic attributes.  Notably, the ggplot theme_grey and theme-bw can be used
#' to change base_size but annoyingly reset my custom legend locations. So apply
#' baseFont last and use it to change just the font sizes without otherwise
#' altering you plot theme.
#'
#' Optionally, set the font family as well. But from what I can tell,
#' on Windows, you can have any font you want as long as it's Helvetica.
#' see: \url{https://github.com/wch/extrafont}
#' regarding attempting to use other font families.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @param base_size Size for the basefont in points (Default = 12)
#' @param base_family Set the font family
#'
#' @return A ggplot theme that can be added to a plot object
#'
#' @examples
#' Myggplot <- myggplot + baseFont(18)
#'
#' @import ggplot2
#'
#' @export
baseFont <- function(base_size=12, base_family=""){
  baseFont <- theme(text = element_text(size=base_size, family = base_family))
}


#some shortcuts for common tweaks I spend alot of time looking up each time.

### Function NoXlab ###
#' Function  NoXlab
#'
#' A simple function to spare me from looking up the syntax everytime
#' I want to remove the X axis label.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @return A ggplot theme that can be added to a plot object
#'
#' @examples
#' Myggplot = myggplot + NoXlab()
#'
#' @import ggplot2
#'
#' @export
NoXlab <- function(){
	NoXlab <- theme(axis.text.x = element_blank())
}

### Function NoYlab ###
#' Function  NoYlab
#'
#' A simple function to spare me from looking up the syntax everytime
#' I want to remove the Y axis label.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @return A ggplot theme that can be added to a plot object
#'
#' @examples
#' Myggplot = myggplot + NoYlab()
#'
#' @import ggplot2
#'
#' @export
NoYlab <- function(){
	NoYlab <- theme(axis.text.y = element_blank())
}

### Function NoLegend ###
#' Function  NoLegend
#'
#' A simple function to spare me from looking up the syntax everytime
#' I want to remove the legend.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @return A ggplot theme that can be added to a plot object
#'
#' @examples
#' Myggplot = myggplot + NoLegend()
#'
#' @import ggplot2
#'
#' @export
NoLegend <- function(){
	NoLegend <- theme(legend.position="none")
}

### Function NoXGrid ###
#' Function  NoXGrid
#'
#' A simple function to spare me from looking up the syntax everytime
#' I want to remove the X grid.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @return A ggplot theme that can be added to a plot object
#'
#' @examples
#' Myggplot = myggplot + NoXGrid()
#'
#' @import ggplot2
#'
#' @export
NoXGrid <- function(){
  NoXGrid <- theme (
  	panel.grid.minor.x=element_blank(),
		panel.grid.major.x=element_blank()
	)
}

### Function NoYGrid ###
#' Function  NoYGrid
#'
#' A simple function to spare me from looking up the syntax everytime
#' I want to remove the Y grid.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @return A ggplot theme that can be added to a plot object
#'
#' @examples
#' Myggplot = myggplot + NoYGrid()
#'
#' @import ggplot2
#'
#' @export
NoYGrid <- function(){
	NoYGrid = theme (
		panel.grid.minor.y=element_blank(),
		panel.grid.major.y=element_blank()
	)
}

### Function setLegendPosition ###
#' Function  setLegendPosition
#'
#' Move the legend around outside (right, left, top, bottom) or inside corner
#' (ne, se, nw, sw) of your plot.  Respects the bw or grey theme style.  You
#' could also use this to get a grey legend on a bw theme or vice versa.
#'
#' This should be called after setting one of greyTheme, bwTheme,
#' theme_grey or theme_bw as those theme reset the legend attributes.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @param ggplot A ggplot plot object
#' @param legendPosition One of top, bottom, right, left, ne, se, nw, sw.  The
#'   first four place it outside the figure. The last four place it inside the
#'   figure. (Default = "right")
#' @param themeStyle one of "grey" or "bw". (Default = "grey")
#'
#' @return A ggplot object with the legend adjusted as specified
#'
#' @examples
#' Myggplot = myggplot + NoXGrid()
#'
#' @import ggplot2
#'
#' @export
setLegendPosition <- function(ggplot, legendPosition = "right", themeStyle = "grey")
{

  # 8 legend positions, with/without grey background
  if (tolower(legendPosition) %in% c("top", "bottom", "left", "right")) {

      ggplot = ggplot +
        theme(legend.position=tolower(legendPosition))

  } else { #place legend inside figure; with or without a grey background

       #four corners
      if (tolower(legendPosition) == "ne") {
        ggplot = ggplot +
          theme(legend.justification=c(1,1), legend.position=c(1,1))
      } else if (tolower(legendPosition) == "se"){
        ggplot = ggplot +
          theme(legend.justification=c(1,0), legend.position=c(1,0))
      } else if (tolower(legendPosition) == "sw"){
        ggplot = ggplot +
          theme(legend.justification=c(0,0), legend.position=c(0,0))
      } else if (tolower(legendPosition) == "nw"){
        ggplot = ggplot +
          theme(legend.justification=c(0,1), legend.position=c(0,1))
      }

  }

  #Set Legend background
  if (tolower(themeStyle) == "bw") {
    ggplot = ggplot +
      theme(legend.background = element_rect(color="grey30", fill="white"))
  } else { #Grey style
    ggplot = ggplot +
      theme(legend.background = element_rect(color="grey30", fill="grey90"))
  }

  return(ggplot)

}

### Function getXrange ###
#' Function  getXrange
#'
#' A simple function to spare me from looking up the syntax everytime
#' I want to get the x range of a plot
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2
#'
#' @param p A ggplot object
#'
#' @return a vector of length 2 with xmin, xmax
#'
#' @examples
#' xrange<- getXrange(Myggplot)
#'
#' @import ggplot2
#'
#' @export
getXrange <- function(p){
  #http://stackoverflow.com/questions/7705345/how-can-i-extract-plot-axes-ranges-for-a-ggplot2-object
  xrange <- ggplot_build(p)$layout$panel_ranges[[1]]$x.range
  #xrange <- ggplot_build(p)$panel$ranges[[1]]$x.range
  xrange <- unlist(xrange)
  return(xrange)
}

### Function getYrange ###
#' Function  getYrange
#'
#' A simple function to spare me from looking up the syntax everytime
#' I want to get the y range of a plot
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2
#'
#' @param p A ggplot object
#'
#' @return a vector of length 2 with xmin, xmax
#'
#' @examples
#' xrange<- getXrange(Myggplot)
#'
#' @import ggplot2
#'
#' @export
getYrange <- function(p){
  #http://stackoverflow.com/questions/7705345/how-can-i-extract-plot-axes-ranges-for-a-ggplot2-object
  yrange <- ggplot_build(p)$layout$panel_ranges[[1]]$y.range
  #yrange <- ggplot_build(p)$panel$ranges[[1]]$y.range
  yrange <- unlist(yrange)
  return(yrange)
}
