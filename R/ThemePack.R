#JRT 1Dec2015
#
#Theme Pack Overview
#
#greyTheme and bwTheme are intended as replacements for the default theme_grey
#and theme_bw themes and are mostly clones of those themes with some modest
#changes. Hadley Wickham defined theme_grey() and theme_bw() functions to
#establish a set of base formatting for plots. theme_grey() is the default
#ggplot theme. theme_bw() just removes the standard grey background but leaves
#other colors alone.
#
#I didn't like how legends scaled with reletive font size changes so implemented
#a legend scale factor that is on by default but can be turned off if you like.
#This works especially well, I think, when the legend is located inside the
#figure borders.
#
#I personally like a little lighter font for the number scales on the X and Y
#axes. But I think Hadley's choice is too light and darkened it a bit (to grey33).
#
#So I cloned theme_grey() and theme(bw) as greyTheme and bwTheme and certain
#attributes to accomplish the above features.
#
#For most purposes you should use greyTheme or bwTheme as your last layer. Since
#these themes reset legend location, you need to change the legend location
#after using one of these thems.  Finally, you can use baseFont to change just
#the base font size without disturbing other features of you plot in order to
#get a font size appropriate for the medium.  I suggest 12 point for knitr
#reports and 24 point for powerpoint presentations.
#
#baseTheme() is a basic theme that sets relative fonts sizes as follows:
#  axis scale numbering = 1.0
#  axis titles = 1.25
#  plot title = 1.5
#  legend text = 1.0
#  legend title = 1.2
#
#Since I incorporated the same relative sizes in greyTheme and bwTheme. There
#is little reason to use this theme. But it does have the advantage that it will
#apply the relative fonts attributes without modifying other customizations.
#
#The function printAndSave, takes the name of a
#plot object and a filename and prints it to the console with small
#and saves it to a file with a larger font. For printAndSave to work properly,
#you should use greyTheme, bwTheme or baseTheme to establish relative font sizes.
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


### Function greyTheme ###
#' Function  greyTheme
#'
#' A customization of the theme_grey default ggplot theme.  Use this instead of
#' theme_grey as the last layer on your plots.  This theme adds relative sizes
#' for most elements (title, axis labels and legends) so that it is easy to
#' adjust all font sizes without otherwise altering your plot. After applying
#' the greyTheme, use the baseFont function to alter just font sizes without
#' otherwise altering your plot. I also found that I didn't like the way legends
#' scaled when shifting from small fonts to larger fonts.  The legend tended to
#' grow unpleasingly large with larger fonts, especially when trying to locate
#' the legend inside the figure borders.  So the legend text is scaled by this
#' equation legendfont = 12/base_size.  You can turn this off with the
#' scaleLegend argument. Finally, I also use a darker color for the axis scale
#' labels (grey33 instead of grey50). I liked the idea of having the axis
#' numbering slightly lighter then the axis title but thought HW's choice of
#' grey50 for the axis numbering was too light.
#'
#' Changes from the default theme_grey applied here are noted in source code comments.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @param base_size Size for the basefont in points (Default = 12)
#' @param base_family Set the font family
#' @param scale_legend Scales Legend.Text to 10/base_size lines (Default=TRUE)
#'
#' @return A ggplot theme that can be added to a plot object
#'
#' @examples
#' Myggplot = myggplot + greyTheme(24)
#'
#' @import ggplot2 grid
#'
#' @export
greyTheme <- function (base_size = 12, base_family = "", scale_legend = TRUE)
{
  #This is my attempt to automatically scale legend sizes.  I want 1 line at 12 point
  #and 0.5 lines at 24 point.  So let's try 12/base_size as a scaling factor.
  if (scale_legend == TRUE && base_size > 17) {
    Legend.ScaledSize <- 10/base_size
    Legend.ScaledFont <- 8/base_size
  } else {
    Legend.ScaledSize = 0.9 #theme_grey defaults 1.2, 0.8
    Legend.ScaledFont = 0.7
  }
  half_line <- base_size / 2

  greyTheme <- theme(
    line = element_line(colour = "black", size = 0.5, linetype = 1, lineend = "butt"),
    rect = element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
    text = element_text(family = base_family,
                        face = "plain", colour = "black", size = base_size, hjust = 0.5,
                        vjust = 0.5, angle = 0, lineheight = 0.9,
                        margin = margin(), debug = FALSE),

    axis.line = element_blank(),
    axis.text = element_text(size = rel(0.8), colour = "grey20"),


    #change axis text to a darker grey; add size
    axis.text.x = element_text(margin = margin(t = 0.8 * half_line / 2),
                               vjust = 1, color="grey30", size=rel(1.0)),
    axis.text.y = element_text(margin = margin(r = 0.8 * half_line / 2),
                               hjust = 1, color="grey30", size=rel(1.0)),

    axis.ticks = element_line(colour = "grey20"),
    axis.ticks.length = unit(half_line / 2, "pt"),

    #add size attribute
    axis.title.x = element_text(size=rel(1.25)),
    axis.title.y = element_text(angle = 90, size=rel(1.25),
                                margin = margin(r = 0.8 * half_line, l = 0.8 * half_line / 2)),  #??? line17 in newTheme_grey.R

    #axis.ticks.margin = unit(0.1, "cm"),

    #border around legend (***different from theme_grey)
    legend.background = element_rect(fill="grey95", colour = "grey30", size=0.5, linetype="solid"),

    legend.margin = unit(0.2, "cm"),
    legend.key = element_rect(fill = "grey95", colour = "white"),

    #scaled
    legend.key.size = unit(Legend.ScaledSize, "lines"),

    legend.key.height = NULL,
    legend.key.width = NULL,

    #scaled
    legend.text = element_text(size = rel(Legend.ScaledFont)),

    legend.text.align = NULL,

    #scaled, bold
    legend.title = element_text(size = rel((Legend.ScaledFont*1.2)), face = "bold"),
    #legend.title = element_text(size = rel((Legend.ScaledFont*1.2))),

    #centered
    legend.title.align = 0.5,

    legend.position = "right",
    legend.direction = NULL,
    legend.justification = "center",
    #legend.box = NULL,

    panel.background = element_rect(fill = "grey90", colour = NA),
    panel.border = element_blank(),
    panel.grid.major = element_line(colour = "white"),
    panel.grid.minor = element_line(colour = "grey95", size = 0.25),
    panel.margin = unit(half_line, "pt"),
    panel.margin.x = NULL,
    panel.margin.y = NULL,
    panel.ontop = FALSE,
    strip.background = element_rect(fill = "grey85", colour = NA),
    strip.text = element_text(colour = "grey10", size = rel(0.8)),

    #add rel size
    strip.text.x = element_text(size=rel(1.6),
                                margin = margin(t = half_line, b = half_line)),
    strip.text.y = element_text(size=rel(1.6), angle = -90,
                                margin = margin(l = half_line, r = half_line)),

    strip.switch.pad.grid = unit(0.1, "cm"),
    strip.switch.pad.wrap = unit(0.1, "cm"),
    plot.background = element_rect(colour = "white"),

    #bold, inc size
    #plot.title = element_text(face="bold", size = rel(1.5),
    plot.title = element_text(size = rel(1.5),
                              margin = margin(b = half_line * 1.2)),

    plot.margin = margin(half_line, half_line, half_line, half_line),
    complete = TRUE

  )
  greyTheme <- greyTheme + baseFont(base_size)
}

### Function bwTheme ###
#' Function  bwTheme
#'
#' A customization of the theme_grey default ggplot theme.  This theme establishes
#' relative sizes for most elements (title, axis labels and legends) so that it is
#' easy to adjust all font sizes without otherwise altering you plot by using the
#' baseFont function.  I also use a darker color for the axis scale labels (grey33 instead of grey50).
#' I like the idea of having the axis title slightly darker than the axis numbering but
#' thought HW's choice of grey50 for the axis numbering was too light.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords ggplot2 fontsize
#'
#' @param base_size Size for the basefont in points (Default = 12)
#' @param base_family Set the font family
#' @param scale_legend Scales Legend.Text to 10/base_size lines (Default=TRUE)
#'
#' @return A ggplot theme that can be added to a plot object
#'
#' @examples
#' Myggplot = myggplot + bwTheme(18)
#'
#' @import ggplot2
#'
#' @export
bwTheme <- function (base_size = 12, base_family = "", scale_legend = TRUE)
{
  #This is my attempt to automatically scale legend sizes.  I want 1 line at 12 point
  #and 0.5 lines at 24 point.  So let's try 12/base_size as a scaling factor.
  if (scale_legend == TRUE && base_size >17) {
    Legend.ScaledSize <- 14/base_size
    Legend.ScaledFont <- 12/base_size
  } else {
    Legend.ScaledSize = 1.2 #theme_grey defaults
    Legend.ScaledFont = 0.8
  }

  greyTheme(base_size = base_size, base_family = base_family) %+replace%
    theme(
      axis.text = element_text(size = rel(0.8)),
      axis.ticks = element_line(colour = "black"),
      legend.key = element_rect(colour = "grey80"),
      panel.background = element_rect(fill = "white", colour = NA),
      panel.border = element_rect(fill = NA, colour = "grey50"),
      panel.grid.major = element_line(colour = "grey90", size = 0.2),
      panel.grid.minor = element_line(colour = "grey98", size = 0.5),
      strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2)
    )
}


### Function baseTheme ###
#' Function  baseTheme
#'
#' A basic theme for individual plots that only sets relative font sizes
#' for common graphic elements.  Apply baseTheme as your last
#' layer.  In particular, theme_bw, theme_grey, bwTheme or greyTheme
#' will undo some of the changes set by baseTheme.  After applying
#' baseTheme, you can use baseFont to adjust font sizes for different
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


#some shortcuts for common tweak I spend alot of time looking up each time.

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
  xrange <- ggplot_build(p)$panel$ranges[[1]]$x.range
  xrange <- unlist(xrange)
  return(xrange)
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
getYrange <- function(p){
  #http://stackoverflow.com/questions/7705345/how-can-i-extract-plot-axes-ranges-for-a-ggplot2-object
  yrange <- ggplot_build(p)$panel$ranges[[1]]$y.range
  yrange <- unlist(yrange)
  return(yrange)
}
