### Function calcPower ###
#' Function  calcPower
#'
#' calcPower is a wrapper around the rnapower function (RNASeqPower package).  
#' calcPower returns a datastructure with power calculations over a range of
#' effect sizes, foldchange and significance levels (alpha).  It takes the
#' coefficient of variation from the matrix of data provided for the analysis.
#' See plotPower for companion functions to plot the data various ways.
#'
#' Provide a numeric data.frame or matrix with two groups of data (e.g. baseline and numerator groups for a ratio)
#'
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords RNA-Seq; Power
#'
#' @param dat A data.frame or matrix of numeric values (CPM, TPM, FPKM or other (non-logged) matrix of data you wish to evaluate).  Columns represent samples and rows represent observations.
#' @param block A vector the same length as ncol(Dat) which identifies which columns represent replicates of a group.
#' @param effect A vector of effect sizes you wish to evaluate (default = seq(1,1, 5.0, 0.1))
#' @param alpha A vector of significance levels you wish to evaluate (default = (0.05, 0.9, 0.05))
#' @param reps A vector of the number of replicates you wish to evaluate (default = c(3, 6, 10))
#'
#' @return A dataframe of power calculation results
#'
#' @examples
#' #run defaults
#' MyContsratList  = runContrasts (MySLOA, ConstrastList)
#' MyContsratList  = runContrasts (MySLOA, ConstrastList, runTopTable = TRUE
#'        ConstrastType = "TopTreat", FoldChangeThreshold = 1.25)
#'
#' @import S4Vectors SummarizedExperiment zFPKM dplyr limma gridExtra
#'
#' @export
calcPower <- function(dat, block, effect = seq(1,1, 5.0, 0.1),
							 alpha=seq(0.05, 0.9, 0.05), 
							 depth=c(10, 100, 1000),
							 reps=c(3, 6, 10){
#set params for power calc
#Liver Power Calculations N=10
#effect <- seq(1.1, 5.0, 0.1)  #c(1.2, 1.5, 2.0, 3.0)
#alpha <- seq(0.05, 0.9, 0.05)  #c(0.05, 0.1, 0.2, 0.5)
#depth <- c(10, 100, 1000)  #1000 reads is approx a median level of expression
#reps=c(3, 6, 10) #number of replicates

#reality checks
	#block must be two groups
	#dat must be all numeric
								 
	#calculate group stats on cpm
	groups <- unique(block)
	groupstats <- list()
	
	#calculate summary stats for each group
	for (group in groups){
	  thisgroup <- list()
	  groupidx <- block == group
	  groupcpm <- dat[,groupidx]
	  #groupcounts <- assay(pRSE, "Counts")[,groupidx]
	  
	  thisgroup[["meanCPM"]] <- rowMeans(groupcpm)
	  #thisgroup[["meanCounts"]] <- rowMeans(groupcounts)
	  thisgroup[["sd"]] <- rowSds(groupcpm)
	  thisgroup[["cv"]] <- thisgroup$sd / thisgroup$meanCPM
	  thisgroup[["meanCoverage"]] <- thisgroup$meanCounts * readlength / exonlength
	  groupstats[[group]] <- thisgroup
	  assign(group, as.data.frame(do.call(cbind, groupstats[[group]])))
	}
	
	
	cv1 = median
	
	cv=median(liver_vehicle$cv[!is.na(liver_vehicle$cv)])
	cv2=median(`liver_BMT-204950`$cv[!is.na(`liver_BMT-204950`$cv)])
	
	
	#core power calculations code
	plist <- list()
	
	for (d in depth){
	  for (n in reps) {
	    rpower <- rnapower (d, n=n, cv=cv, cv2=cv2, effect=effect, alpha = alpha) %>% 
	      as.data.frame %>% round(3)
	    #colnames(rpower) = paste("FDR",colnames(rpower), sep="_")
	    rpower$EffectSize <- rownames(rpower)
	    rpower$N <- n
	    rpower$Depth <- d
	    dfname <- paste("power_", "D", d, "_N", n, sep="")
	    plist[[dfname]] <- rpower #assign(dfname, rpower)
	    # assign (dfname, rpower)
	    # plist[[dfname]] <- rpower
	  }
	}
	
	#concatenate the power tables
	allpower <- (do.call("rbind", plist) )
	tallpower <- gather(allpower, fdr, power, -c(EffectSize, N, Depth))
	tallpower$EffectSize %<>% as.numeric
	tallpower$fdr %<>% as.numeric
	
	
	#Plots
	
	#ROC curves
	
	
	#N=10 data
	rocdat <- filter (tallpower, N==10, EffectSize %in% c(1.2, 1.5, 2, 3))
	rocdat$Depth %<>% as.factor
	rocdat$EffectSize %<>% as.factor
	rocdat$EffectSize <- paste("FoldChg", rocdat$EffectSize,sep=" = ")
	rocdat$N %<>% as.factor
	
	roc <- ggplot(rocdat, aes(x=fdr, y=power, fill=Depth, shape=Depth, color = Depth)) +
	  geom_point(size = 2) +
	  geom_line(size = 1) +
	  scale_y_continuous(breaks=seq(0,1,0.2), minor_breaks=seq(0,1,0.05)) +
	  scale_x_continuous(breaks=seq(0,1,0.2)) +
	  facet_wrap(~ EffectSize, ncol=2) +
	  geom_vline(xintercept=0.1, color="red") +
	  ggtitle("ROC Curves for Mouse Liver Data (n=10)") +
	  xlab ("\nFDR") +
	  ylab ("Power") +
	  greyTheme()
	
	printAndSave(roc, "PowerROC.PNG")
	
	#Power as a function of EffectSize (FoldChange)
	
	#N=10 data, fdr 10%
	rocdat <- filter (tallpower, N==10, fdr == 0.1) 
	rocdat$Depth %<>% as.factor
	# rocdat$EffectSize %<>% as.factor
	# rocdat$EffectSize <- paste("FoldChg", rocdat$EffectSize,sep=" = ")
	# rocdat$N %<>% as.factor
	
	roc <- ggplot(rocdat, aes(x=EffectSize, y=power, fill=Depth, shape=Depth, color = Depth)) +
	  geom_point(size = 3) +
	  geom_line(size = 1) +
	  scale_x_continuous(breaks=seq(1,5,0.5)) +
	  scale_y_continuous(breaks=seq(0,1,0.1)) +
	 #facet_wrap(~ EffectSize, ncol=2) +
	  geom_vline(xintercept=c(1.5,2), color="red") +
	  ggtitle("Power Curves at N=10, 10% FDR") +
	  xlab("\nFold Change") +
	  ylab("Power") +
	  greyTheme()
	
	printAndSave(roc, "PowerVsEffectSize.PNG")
	
	#Power as a function of Replicates
	
	#N=10 data, fdr 10%
	rocdat <- filter (tallpower, N==10, fdr == 0.1) 
	rocdat$Depth %<>% as.factor
	# rocdat$EffectSize %<>% as.factor
	# rocdat$EffectSize <- paste("FoldChg", rocdat$EffectSize,sep=" = ")
	# rocdat$N %<>% as.factor
	
	roc <- ggplot(rocdat, aes(x=EffectSize, y=power, fill=Depth, shape=Depth, color = Depth)) +
	  geom_point(size = 3) +
	  geom_line(size = 1) +
	  scale_x_continuous(breaks=seq(1,5,0.5)) +
	  scale_y_continuous(breaks=seq(0,1,0.1)) +
	 #facet_wrap(~ EffectSize, ncol=2) +
	  geom_vline(xintercept=c(1.5,2), color="red") +
	  ggtitle("Power Curves at N=10, 10% FDR") +
	  xlab("\nFold Change") +
	  ylab("Power") +
	  greyTheme()
	
	printAndSave(roc, "PowerVsEffectSize.PNG")
	
	

