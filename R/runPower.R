### Function runPower ###
#' Function runPower (DGE.Tools2)
#'
#' Take a counts matrix and design matrix and return a power analysis using
#' the RNASeqPower package.  The counts matrix should be prefiltered to remove
#' non-expressed genes by your favorite filtering critera.  The design matrix
#' should describe the major sources of variation so the procedure can dial
#' out those known effects for the power calculations.
#'
#' If return = "dataframe" the function will return a tall skinny dataframe
#' of power calculation for various requested combinations of N and signficance
#' thresholds.  If return = "plots" or "both", a list is returned with two
#' ggplots (plots) or the plots plus the dataframe (both).
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq, DGEobj
#'
#' @param counts  A counts matrix (required)
#' @param designMatrix A design matrix (required)
#' @param depth A  set of depth to use in the calculations.  The default depths of
#' c(10, 100, 1000) respectively represent a detection limit, below average
#' expression and median expression levels, express in readcount units.
#' @param N A set of N value to report power for (default = c(3, 6, 10, 20))
#' @param effectSize A set of fold change values to test (default = c(1.2, 1.5, 2))
#' @param return One of "dataframe", "plots", "both" (default = "both").
#' Two plots are generated; a ROC curve (FDR vs Power) and a plot of N vs Power.
#'
#' @return A DGEobj with a modified object definition embedded.
#'
#' @examples
#'    MyResults <- runPower(counts, designMatrix)
#'
#' @import magrittr assertthat RNASeqPower
#'
#' @export
runPower <- function(counts, designMatrix,
                     depth = c(10, 100, 1000),
                     N = c(3, 6, 10, 20),
                     FDR = c(0.05, 0.1),
                     effectSize = c(1.2, 1.5, 2),
                     return = "both"){

#Fit the BCV data and define the BCV for each depth requested.
  #estimate dispersion
  dgelist <- counts %>%
    as.matrix %>%
    DGEList %>%
    calcNormFactors %>%
    estimateDisp (design=designMatrix, robust=TRUE)

  # BCV <- sqrt(dgelist$tagwise.dispersion)

  ### Get a fitted CV values for each input value of depth
  #need to calculate the cv given a depth value
  #estimateDisp add trended.dispersion to the dgelist.
  #But this was calculated against the AveLogCPM.
  #Need to convert the nominal counts (define by the depth argument)
  # to AveLogCPM units and return the BCV for those values
  #BCV is the sqrt of Dispersion
  GeoMeanLibSize <-  dgelist$samples$lib.size %>% log %>% mean %>% exp
  depth_avelogcpm <- aveLogCPM(depth, GeoMeanLibSize)
  depthBCV <- sqrt(approx(dgelist$AveLogCPM, dgelist$trended.dispersion,
                               xout=depth_avelogcpm, rule=2)$y)


  #generate a tall skinny dataframe with all combinations of depth, n, effect, alpha and calculated power


  #set other values

  n <- seq(min(N),max(N),1)   #for the N vs P plot
  alpha <- seq(0.05, 0.9, 0.05)  #alpha is FDR levels

  #initialize a dataframe for the results table
  pdat <- data.frame(depth=double(),
                     n = double(),
                     effect = double(),
                     alpha = double(),
                     power = double(),
                     stringsAsFactors=FALSE)
  cnames <- colnames(pdat)

  for (D in depth){
    cv <- depthBCV[D==depth]
    for (Nf in n)
      for (E in effectSize)
        for (A in alpha) {
          P <- rnapower(depth=D, n=Nf, cv=cv, effect=E, alpha=A)
          pdat <- rbind(pdat, c(depth=D, n=Nf, effect=E, alpha=A, power=P))
        }
  }
  colnames(pdat) <- cnames



  result <- list()
  if (tolower(return) %in% c("dataframe", "both")){
    result$PowerData <- pdat
  }

  if (tolower(return) %in% c("plot", "both")){
    #draw the plots
    #
    #ROC Curves (FDR vs Power)

    rocdat <- filter(pdat, n %in% N)
    rocdat$depth %<>% as.factor

    roc <- ggplot(rocdat, aes(x=alpha, y=power, fill=depth, shape=depth, color = depth)) +
      # geom_point(size = 3) +
      geom_line(size = 1) +
      scale_x_continuous(breaks=seq(0, 1, 0.2)) +
      scale_y_continuous(breaks=seq(0,1,0.2)) +
      facet_grid(effect ~ n, labeller = label_both) +
      # geom_vline(xintercept=0.1, color="red") +
      ggtitle("ROC curves") +
      xlab("\nFDR") +
      ylab("Power") +
      expand_limits(x=0,y=0) +
      theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      theme_gray(18)

    result$ROC <- roc

    #N vs Power
    #filter to just a few FDR thresholds
    ndat <- filter(pdat, alpha %in% FDR)

    ndat$depth %<>% as.factor
    ndat$FDR <- ndat$alpha
    # ndat$effect %<>% as.factor
    # ndat$n %<>% as.factor

    NvP <- ggplot(ndat, aes(x=n, y=power, fill=depth, shape=depth, color = depth)) +
      # geom_point(size = 3) +
      geom_line(size = 1) +
      # scale_x_continuous(breaks=seq(0, 1, 0.2)) +
      scale_y_continuous(breaks=seq(0,1,0.2)) +
      facet_grid(FDR ~ effect, labeller = label_both) +
      ggtitle("N vs Power") +
      xlab("\nN") +
      ylab("Power") +
      expand_limits(x=0,y=0) +
      theme_gray()
    
    result$NvP <- NvP
  }
  return(result)
}




