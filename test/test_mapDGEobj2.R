library(igraph)
library(RColorBrewer)

# Prepare an iGraph object for plotting
mynet <- mapDGEobj(dgeObj)

#define a layout
lay.sug <- layout_with_sugiyama(mynet)
lay.tree <- layout_as_tree

#2D Plot
plot(mynet,
     edge.arrow.size=.3,
     vertex.color="dodgerblue2",
     vertex.frame.color="dodgerblue4",
     layout=lay.sug$layout)

#2D Plot; colorby the basetype attribute
pal <- brewer.pal(length(unique(V(mynet)$basetype)), "Set1")
myPallet <- pal[as.numeric(as.factor(vertex_attr(mynet, "basetype")))]

plot(mynet,
     edge.arrow.size=.3,
     vertex.color = myPallet,
     vertex.label.family = "Helvetica",
     layout=lay.sug$layout,
     edge.arrow.size = 3,
     edge.arrow.width = 3,
     edge.color = "black",
     edge.size = 3,
     vertex.size = 20,
     # vertex.shape = "rectangle",
     vertex.size2 = 20,
     vertex.label.dist = -pi/2
)

#2D Interactive plot
plotHandle <- tkplot(mynet,
                     vertex.color=myPallet,
                     vertex.frame.color="dodgerblue4",
                     canvas.width = 800,
                     canvas.height = 800,
                     layout=lay.sug$layout)
tk_close(plotHandle)
