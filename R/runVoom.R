### Function runVoom ###
#' Function  runVoom (voom then lmFit)
#'
#' In the recommended workflow this runs voomWithQualityWeights followed by
#' lmfit and optionally eBayes.  You should enable eBayes if the contrasts of
#' interest are already represented in the model. If you intend to use
#' contrasts.fit downstream, you should run eBayes after that step instead. In
#' other words, run eBayes last.
#'
#' Input is minimally a DGEobj containing a DGEList (typically TMM-normalized),
#' and a formula (text representation).  Other arguments can invoke
#' duplicateCorrelation and modify use of quality weights.
#'
#' Returns a DGEobj class object containing the designMatrix, VoomElist (voom
#' output) and Fit object (lmfit output). Appends data items to the input
#' DGEobj.
#'
#' Quality weights should be left enabled unless you have a good reason to turn it
#' off.  If all samples are equal quality, the weights will all approach 1.0 with no
#' consequence on the results.  More typically,  some samples are better than others
#' and using quality weights improves the overall result.
#'
#' Use var.design when you notice that quality weights are correlated with some
#' factor in the experiment.  This will cause the quality weights to be
#' calculated as a group instead of individually.
#'
#' Use duplicate correlattion when you have subjects that have been sampled more
#' than once (e.g. before and after some treatment).  This calculates a within
#' subject correlation and includes this in the model.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID
#'
#' @param dgeObj A DGEobj containing a DGEList (e.g. from runEdgeRNorm) (required)
#' @param designMatrixName Name of a design matrix within dgeObj (required)
#' @param dupcorBlock Supply a block argument to trigger duplicateCorrelation (optional).
#'    Should be a vector the same length as ncol with values to indicate common
#'    group membership for duplicateCorrelation.
#' @param qualityWeights Runs VoomWithQualityWeights if set to TRUE (default=TRUE).
#'    This should normally be set to TRUE.
#' @param var.design Provide a design matrix (from model.matrix) to identify
#'    replicate groups (e.g. "~ ReplicateGroup") for quality weight determination.
#'    Causes quality weights to be determined on a group basis.  If omitted
#'    VoomWithQualityWeights treats each sample individually.
#' @param runEBayes Runs eBayes after lmFit; default = TRUE
#' @param robust Used by eBayes (Default = TRUE)
#' @param proportion Proportion of genes expected to be differentially expressed
#'   (used by eBayes) (Default = 0.01) Modify the prior accordingly if your 1st pass analysis shows
#'   a significantly higher or lower proportion of genes regulated than the default.
#' @param mvPlot Enables the voom mean-variance plot (Default = TRUE)
#' @return A DGEobj now containing designMatrix, Elist and fit object
#'
#' @examples
#'
#' @import magrittr
#' @importFrom limma voom lmFit eBayes voomWithQualityWeights duplicateCorrelation
#' @importFrom stringr str_c
#' @importFrom DGEobj getItem addItem
#' @importFrom assertthat assert_that
#'
#' @export
runVoom <- function(dgeObj, designMatrixName,
                    dupcorBlock,
                    qualityWeights = TRUE,
                    var.design,
                    mvPlot = TRUE,
                    runEBayes = TRUE,
                    robust = TRUE,
                    proportion=0.01
                    ){

    assert_that(!missing("dgeObj"),
                !missing("designMatrixName"),
                designMatrixName %in% names(dgeObj),
                class(dgeObj)[[1]] == "DGEobj")

    designMatrix <- DGEobj::getItem(dgeObj, designMatrixName)

    #get the DGEList
    if ("DGEList" %in% attr(dgeObj, "type"))
        dgelist <- DGEobj::getItem(dgeObj, "DGEList")
    else stop("No DGEList found in DGEobj")

    #collect calling args for documentation
    funArgs <- match.call()

    #Set run parameters
    dupcor <- FALSE
    if (!missing(dupcorBlock))
        dupcor <- TRUE

    blockQW <- FALSE
    if (qualityWeights == TRUE & !missing(var.design))
        blockQW <- TRUE

    #option flags  dupcor, qualityWeights, blockQW
    # FFF: "baseOptions"
    # FTF  "indQW"
    # FTT  "blockedQW"
    # TFF  "dupcor_base"
    # TTF  "dupcor_indQW"
    # TTT  "dupcor_vdQW" (var.design blocking on QW)

# Note re Duplicate Correlation:
# Gordon Smythe has noted that duplicateCorrelation should be run twice.
# Actually, should interate until asymtope achieved but normally twice is sufficient.
# https://support.bioconductor.org/p/59700/#67620

#Main Calculations (one of six blocks will be run)

    #set type of analysis
    if (dupcor==F & qualityWeights==F & blockQW==F){ #baseOptions

        #voom squeezes the variance (borrowing from other genes) to deal
        #with the heteroskedasticity problem
        VoomElist <- limma::voom(dgelist, designMatrix, plot=mvPlot, col="blue")
        fit <- limma::lmFit(VoomElist, designMatrix)

    } else if (dupcor==F & qualityWeights==T & blockQW==F){ #indQW analysis

        VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix,
                                            plot=mvPlot, col="blue")
        fit <- limma::lmFit(VoomElist, designMatrix)

    } else if (dupcor==F & qualityWeights==T & blockQW==T){ #blockedQW analysis

        VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix,
                               plot=mvPlot, col="blue",
                               var.design = var.design)
        fit <- limma::lmFit(VoomElist, designMatrix)

    } else if (dupcor==T & qualityWeights==F & blockQW==F){ #dupcor_base analysis

        VoomElist <- limma::voom(dgelist, designMatrix)
        corfit <- limma::duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)
        VoomElist <- limma::voom(dgelist, designMatrix,
                          correlation=corfit$consensus.correlation,
                          plot=mvPlot, col="blue")
        corfit <- limma::duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)

        fit <- limma::lmFit(VoomElist, designMatrix, block=dupcorBlock,
                  correlation=corfit$consensus.correlation)

    } else if (dupcor==T & qualityWeights==T & blockQW==F){ #dupcor_indQW analysis

        VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix)
        corfit <- limma::duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)
        VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix,
                                            plot=mvPlot, col="blue",
                                            correlation=corfit$consensus.correlation)
        corfit <- limma::duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)

        fit <- limma::lmFit(VoomElist, designMatrix, block=dupcorBlock,
                     correlation=corfit$consensus.correlation)

    } else if (dupcor==T & qualityWeights==T & blockQW==T){ #dupcor_vdQW analysis

        VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix,
                                            var.design = var.design)
        corfit <- limma::duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)
        VoomElist <- limma::voomWithQualityWeights(dgelist, designMatrix,
                                            plot=mvPlot, col="blue",
                                            correlation=corfit$consensus.correlation,
                                            var.design = var.design)
        corfit <- limma::duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)

        fit <- limma::lmFit(VoomElist, designMatrix, block=dupcorBlock,
                     correlation=corfit$consensus.correlation)

    } else {
        stop ("argument combination not supported!")
    }

    #run eBayes
    if (runEBayes) {
        fit = limma::eBayes(fit, robust=robust, proportion=proportion)
        itemAttr <- list(eBayes = TRUE)
    } else itemAttr <- list(eBayes = FALSE)

    if (exists("corfit")) { #duplicate correlation was used; capture the correlation value
      cat (stringr::str_c("Duplicate Correlation = ", round(corfit$consensus.correlation, 4), "   \n"))
      attr(VoomElist, "DupCor") <- corfit$consensus.correlation
      attr(fit, "DupCor") <- corfit$consensus.correlation
    }

    #save the several objects

    VoomElistName = paste(designMatrixName, "_Elist", sep="")
    dgeObj %<>% DGEobj::addItem(VoomElist, VoomElistName,
                        "Elist",
                        funArgs=funArgs,
                        parent=list("DGEList", designMatrixName)
                        )

    #Add corfit if present
    if (exists("corfit"))
      dgeObj %<>% DGEobj::addItem(corfit, paste(designMatrixName, "_corFit", sep=""),
                                  "corFit",
                                  funArgs=funArgs,
                                  parent=paste(designMatrixName, "_Elist", sep=""))

    dgeObj %<>% DGEobj::addItem(fit, paste(designMatrixName, "_fit", sep=""),
                                "fit",
                                funArgs=funArgs,
                                itemAttr=itemAttr,
                                parent=list(VoomElistName, designMatrixName))

  return(dgeObj)  #now containing designMatrix, corFit, Elist and fit
}
