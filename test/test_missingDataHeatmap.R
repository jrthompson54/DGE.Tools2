#test missingDataHeatmap

#make up some data with one row of NAs
df <- data.frame(Doubles=c(1:10, NA, 12:20),
                       Characters=c(rep("A",10), NA, rep("B",9)),
                       stringsAsFactors=FALSE)
rownames(df) <- paste("A", rownames(df), sep="")

missingDataHeatmap(df, title="Test Title", legend=FALSE, sub="Red indicates missing data")

missingDataHeatmap(df, title="test", fontsize=10, color="RdYlBu2")


