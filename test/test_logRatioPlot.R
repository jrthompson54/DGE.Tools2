x <- structure(list(EnsgID = c("Gene1", "Gene2", "Gene3", "Gene4",
                          "Gene1", "Gene2", "Gene3", "Gene4", "Gene1", "Gene2", "Gene3",
                          "Gene4", "Gene1", "Gene2", "Gene3", "Gene4"),
               logFC = c(4.245, 2.428, 1.892, 6.746, -0.944, 0.16, 0.131, -1.223, -1.96, -0.973,
                         0.093, -3.144, -1.607, 0.068, 0.216, -1.614),
               CI.L = c(3.558, 2.01, 1.452, 6.08, -1.173, -0.006, -0.075, -1.495, -2.416, -1.277,
                        -0.213, -3.781, -1.85, -0.093, 0.024, -1.888),
               CI.R = c(4.933, 2.846, 2.333, 7.413, -0.716, 0.326, 0.337, -0.951, -1.504, -0.669,
                        0.4, -2.507, -1.364, 0.229, 0.408, -1.341),
               Contrast = c("Contrast1", "Contrast1", "Contrast1", "Contrast1", "Contrast2", "Contrast2",
                            "Contrast2", "Contrast2", "Contrast3", "Contrast3", "Contrast3",
                            "Contrast3", "Contrast4", "Contrast4", "Contrast4", "Contrast4"
               )), row.names = c(NA, -16L), class = "data.frame")
#facetted barplot
MyPlot <- logRatioPlot(x,
                       facetColname = "EnsgID",
                       xColname = "Contrast",
                       facetCol = 4,
                       scales = "fixed",
                       barWidth = 0.7)
assertthat::assert_that("ggplot" %in% class(MyPlot))

#facetted point plot
MyPlot2 <- logRatioPlot(x, plotType = "point",
                       facetColname = "EnsgID",
                       xColname = "Contrast",
                       facetCol = 4,
                       scales = "fixed")
assertthat::assert_that("ggplot" %in% class(MyPlot2))

#individual point plots
MyPlot3 <- logRatioPlot(x, plotType = "point",
                        facetColname = "EnsgID",
                        xColname = "Contrast",
                        facet = FALSE,
                        scales = "fixed")
assertthat::assert_that(class(MyPlot3 == "list"))
