library(tidyverse)
library(DGE.Tools2)

#create a dataframe of data with row and colnames
groups <- paste("group", factor(rep(1:4, each=100)), sep="")
x <- matrix(rnorm(2400, mean=10), ncol=length(groups))
colnames(x) <- paste("sample", 1:ncol(x), sep="")
rownames(x) <- paste("gene", 1:nrow(x), sep="")
#reformat into tidy dataframe
tidyInt <- tidyIntensity(x,
                         rowidColname ="GeneID",
                         keyColname = "Sample",
                         valueColname ="Log2CPM",
                         group=groups)
#Facetted boxplot
obsPlot2(tidyInt, plotByCol="GeneID",
         groupCol="group",
         valueCol ="Log2CPM",
         pointJitter = 0.1,
         facetRow= 2)
#Facetted violin plot
obsPlot2(tidyInt, plotByCol="GeneID",
         violinLayer = TRUE,
         boxLayer = FALSE,
         groupCol="group",
         valueCol ="Log2CPM",
         pointJitter = 0.1,
         facetRow= 2)
