name
#try all 6 senarios
d1 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = FALSE)

#QW and Var.design
vd <- model.matrix(as.formula("~ Treatment"), design)
d2 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = FALSE,
              var.design=vd)

#QW and Var.design and dupCor
block <- c(1,2,3,1,2,3,4,5,6,4,5,6,7,8,9,7,8,9)
d3 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = FALSE,
              var.design=vd,
              dupcorBlock=block)
#add QW
d4 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = TRUE)

#QW and Var.design
vd <- model.matrix(as.formula("~ Treatment"), design)
d5 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = TRUE,
              var.design=vd)

#QW and Var.design and dupCor
block <- c(1,2,3,1,2,3,4,5,6,4,5,6,7,8,9,7,8,9)
d6 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = TRUE,
              var.design=vd,
              dupcorBlock=block)

#QW and Var.design with eBayes
vd <- model.matrix(as.formula("~ Treatment"), design)
d7 <- runVoom(dgeObj, designMatrix, formula, qualityWeights = TRUE,
              var.design=vd, 
              runEBayes=TRUE)
