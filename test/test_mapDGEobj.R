x <- getwd()
setwd("C:/Users/thompj27/Documents/Fibrosis/IPF_Lung_(UPenn)/P-20161128-0001/RNA-Seq_Analysis/P00931_IPF_Lung_UPenn_2Batches_P-20161128-0001/ClinData")
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)
library(tidyverse)
library(magrittr)
library(igraph)
library(RColorBrewer)
dgeObj <- readRDS("../output/IPF_Lung_UPenn_P-20161128-0001_11Sep2017_DGEobj.RDS")

net <- mapDGEobj(dgeObj)

#apply sugiyama layout style
lay.sug <- layout_with_sugiyama(net)
# setwd("~/R/iGraphTutorial")
#
#
# nodes2 <- read.csv("./Data files/Dataset2-Media-User-Example-NODES.csv", header=T, as.is=T)
# links2 <- read.csv("./Data files/Dataset2-Media-User-Example-EDGES.csv", header=T, row.names=1)


#2dPlot
plot(net,
     edge.arrow.size=.3,
     vertex.color="dodgerblue2",
     vertex.frame.color="dodgerblue4",
     canvas.width = 1000, canvas.height = 2000,
     layout=lay.sug$layout)


#2d interactive plot
tkplot(net, vertex.color="dodgerblue2", vertex.frame.color="dodgerblue4",
       canvas.width = 800, canvas.height = 800,
       layout=lay.sug$layout)


pal <- brewer.pal(length(unique(V(net)$type)), "Set1")
myPallet <- pal[as.numeric(as.factor(vertex_attr(net, "type")))]
plot(net, edge.arrow.size=.3,
     vertex.color = myPallet,
     vertex.label.dist = 2,
     vertex.label.family = "Helvetica")


plotHandle <- tkplot(net, vertex.color=pal[as.numeric(as.factor(vertex_attr(net, "type")))],
       vertex.frame.color="dodgerblue4",
       canvas.width = 800,
       canvas.height = 800)


zach.sug <- layout_with_sugiyama(zach)
  plot(zach,
       edge.arrow.size=.3,
       vertex.color = myPallet,
       vertex.label.family = "Helvetica",
       layout=zach.sug$layout
       )

