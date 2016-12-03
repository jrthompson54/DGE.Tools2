### Function runVoom ###
#' Function  runVoom (voom then lmFit)
#'
#' Input is minimally a DGEList (typically from counts %>% DGEList %>%
#' calcNormFactors), a design matrix (from model.matrix) and a formula (text
#' representation).  Instead of a dgelist, you can probably use a Log2CPM matrix
#' (preferably already normalized).
#' 
#' Returns a dgeResult object containing the VoomElist and Fit objects derived from
#' voom and lmfit respectively. Appends to a prior dgeResult if one is passed, else
#' creates a new dgeResult.
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID
#'
#' @param dgelist An edgeR DGEList object that has already been normalized 
#' (e.g. by edgeR::calcNormFactors or DGE.Tools::runEdgeRNorm) (required)
#' @param designMatrix A design matrix created by model.matrix (required)
#' @param formula A text representation of a formula to fit (required)
#' @param dupcorBlock Supply a block argument to trigger duplicateCorrelation (optional)
#' @param qualityWeights Runs VoomWithQualityWeights if set to TRUE (default=TRUE)
#' @param var.design Provide a design matrix to identify replicate groups.  If omitted
#'    VoomWithQualityWeights treats each sample individually. 
#' @param dgeResult A List object to add results to. If supplied, results 
#'    will be added to this. Otherwise a new dgeResult will be created. This
#'    allows you to sequentially capture results of multiple steps by passing
#'    the dgeResult to each function in the pipeline that supports it. 
#'
#' @return A dgeResult object (currently a simple List of Lists structure)
#'
#' @examples
#'
#' @import magrittr limma 
#'
#' @export
runVoom <- function(dgelist, designMatrix, formula,
                    dupcorBlock=NULL,
                    qualityWeights = TRUE,
                    var.design=NULL,
                    dgeResult=NULL){

  #Version Info
  FVersion = "runVoom : 01Dec2016"

    #check args 
    if (!exists("dgelist"))
        stop("dgelist is a required argument")
    if (!exists("designMatrix"))
      stop("designMatrix is a required argument")
    if (!exists("formula"))
      stop("formula is a required argument")
    
    #check datatypes
    if (!class(designMatrix) == "matrix")
        stop("designMatrix Arg must be type matrix (created by model.matrix)")
    #check formula format
    result <- try(as.formula(formula), silent=TRUE)
    if (class(result) == "try-error") {
      stop("You're forumla is badly formatted")
    } 
    
    #Set run parameters  
    dupcor <- FALSE
    if (!is.null(dupcorBlock))
        dupcor <- TRUE
    
    blockQW <- FALSE
    if (qualityWeights == TRUE & !is.null(var.design))
        blockQW <- TRUE
    
    workflowRecord <- data.frame(Function=as.character, 
                                 Version=as.character,
                                 Date=as.Date(character()),
                                 RunParam=as.character,
                                 stringsAsFactors = FALSE)
    #inherit workflowRecord List from the dgeResult if present
    if(!is.null(dgeResult))
        if(with(dgeResult, exists('workflowRecord'))) 
            workflowRecord <- dgeResult$workflowRecord

    
    #create dgeResult if not provided
    if(is.null(dgeResult))
        dgeResult <- list()
    
    #option flags  dupcor, qualityWeights, blockQW
    # FFF: "baseOptions"
    # FTF  "indQW"
    # FTT  "blockedQW"
    # TFF  "dupcor_base"
    # TTF  "dupcor_indQW"
    # TTT  "dupcor_vdQW" (var.design blocking on QW)

# Note: Gordon Smythe has noted that duplicateCorrelation should be run twice.
# Actually, should interate until asymtope achieved but normally twice is sufficient.

#workflow record: DF with cols FunctionName, Version, Date, RunParams

#Main Calculations (one of six blocks will be run)
    
    #set type of analysis
    if (dupcor==F & qualityWeights==F & blockQW==F){ #baseOptions
        
        #workflow contains: FunName, version, date, run params
        workflow <- c("runVoom", FVersion, date(), 
                                 "dupcor OFF, QW OFF, QWblock OFF")
        #voom squeezes the variance (borrowing from other genes) to deal
        #with the heteroskedasticity problem
        VoomElist <- voom(dgelist, designMatrix)
        #collect VoomElist in results
        dgeResult$VoomElist$data <- VoomElist
        dgeResult$VoomElist$type <- "row"
        fit <- lmFit(VoomElist, designMatrix)
        #collect fit in results
        dgeResult$fit$data <- fit
        dgeResult$fit$type <- "row"
        
    } else if (dupcor==F & qualityWeights==T & blockQW==F){ #indQW analysis
        
        workflow <- c("runVoom", FVersion, date(), 
                      "dupcor OFF, QW ON, QWblock OFF")
        VoomElist <- voomWithQualityWeights(dgelist, designMatrix, plot=FALSE)
        dgeResult$VoomElist$data <- VoomElist
        dgeResult$VoomElist$type <- "row"
        fit <- lmFit(VoomElist, designMatrix)
        
    } else if (dupcor==F & qualityWeights==T & blockQW==T){ #blockedQW analysis
        
        workflow <- c("runVoom", FVersion, date(), 
                      "dupcor OFF, QW ON, QWblock ON")
        VoomElist <- voomWithQualityWeights(dgelist, designMatrix, plot=FALSE,
                               var.design = var.design)
        dgeResult$VoomElist$data <- VoomElist
        dgeResult$VoomElist$type <- "row"
        fit <- lmFit(VoomElist, designMatrix)
        
    } else if (dupcor==T & qualityWeights==F & blockQW==F){ #dupcor_base analysis
        
        workflow <- c("runVoom", FVersion, date(), 
                      "dupcor ON, QW OFF, QWblock OFF")
        VoomElist <- voom(dgelist, designMatrix)
        corfit <- duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)
        VoomElist <- voom(dgelist, designMatrix, correlation=corfit$consensus.correlation)
        corfit <- duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)
        dgeResult$VoomElist$data <- VoomElist
        dgeResult$VoomElist$type <- "row"
        fit <- lmFit(VoomElist, designMatrix, block=dupcorBlock,
                  correlation=corfit$consensus.correlation)
        
    } else if (dupcor==T & qualityWeights==T & blockQW==F){ #dupcor_indQW analysis
        
        workflow <- c("runVoom", FVersion, date(), 
                      "dupcor ON, QW ON, QWblock OFF")
        VoomElist <- voomWithQualityWeights(dgelist, designMatrix, plot=FALSE)
        corfit <- duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)
        VoomElist <- voomWithQualityWeights(dgelist, designMatrix, plot=FALSE,
                                            correlation=corfit$consensus.correlation)
        corfit <- duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)
        dgeResult$VoomElist$data <- VoomElist
        dgeResult$VoomElist$type <- "row"
        fit <- lmFit(VoomElist, designMatrix, block=dupcorBlock,
                     correlation=corfit$consensus.correlation)
        
    } else if (dupcor==T & qualityWeights==T & blockQW==T){ #dupcor_vdQW analysis
        
        workflow <- c("runVoom", FVersion, date(), 
                      "dupcor ON, QW ON, QWblock ON")
        VoomElist <- voomWithQualityWeights(dgelist, designMatrix, plot=FALSE,
                                            var.design = var.design)
        corfit <- duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)
        VoomElist <- voomWithQualityWeights(dgelist, designMatrix, plot=FALSE,
                                            correlation=corfit$consensus.correlation,
                                            var.design = var.design)
        corfit <- duplicateCorrelation(VoomElist,
                                       designMatrix,
                                       block=dupcorBlock)
        dgeResult$VoomElist$data <- VoomElist
        dgeResult$VoomElist$type <- "row"
        fit <- lmFit(VoomElist, designMatrix, block=dupcorBlock,
                     correlation=corfit$consensus.correlation)
        
    } else {
        stop ("argument combination not supported!")
    }

    workflowRecord %<>% rbind(workflow)
    dgeResult$workflowRecord$dat <- workflowRecord
    dgeResult$workflowRecord$type <- "meta"

  return(dgeResult)  #now containing Fit info etc.
}
