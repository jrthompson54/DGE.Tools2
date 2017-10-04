rm(list=ls())
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)



#omicsoft DGEobj
dge_as <- readRDS("../DGEobj.RDS")
counts <- getItem(dge_as, "counts")
geneData <- getItem(dge_as, "geneData")
design <- getItem(dge_as, "design")
tpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize="tmm", debug=T)
logtpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize="tmm", log=TRUE, debug=T)
fpkm <- convertCounts(counts, unit="fpkm", geneLength=geneData$ExonLength, normalize="none", debug=T)
fpkm_tmm  <- convertCounts(counts, unit="fpkm", geneLength=geneData$ExonLength, normalize="tmm", debug=T)
logfpkm <- convertCounts(counts, unit="fpkm", geneLength=geneData$ExonLength, normalize="none", log=TRUE, debug=T)
logfpkm_tmm <- convertCounts(counts, unit="fpkm", geneLength=geneData$ExonLength, normalize="tmm", log=TRUE, debug=T)

inspect(counts)
inspect(tpm)
inspect(logtpm)
inspect(fpkm)
inspect(fpkm_tmm)
inspect(logfpkm)
inspect(logfpkm_tmm)

#what happens to tpm=0 cells (prior should force them to a low log value)
idx <- tpm==0
head(logtpm[idx])

#check distributions (expect medians around 3-5 (logged))
hist(log2(tpm[,1]), breaks=100)
hist(logtpm[,1], breaks=100)
hist(log2(fpkm[,1]), breaks=100)
hist(log2(fpkm_tmm[,1]), breaks=100)
hist(logfpkm[,1], breaks=100)
hist(logtpm[,1], breaks=100)
hist(logfpkm_tmm[,1], breaks=100)

#Show effect of prior.count on logtpm
plot(logtpm[,1], log2(tpm[,1]))
#note the prior.count squeezes the lowest values (-10 to -20) to above -10 otherwise effect is negligible.


#Confirm that TPM is the same as TMM TPM
tmmtpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize="tmm", debug=T)
tpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize="none", debug=T)
plot(log2(tmmtpm[,1]), log2(tpm[,1]))
diff <- tmmtpm - tpm
hist(diff[,1]) #not all zeros but very close.
max(diff)  #biggest diff is 10e-11
inspect(diff)
cor(tmmtpm[,1], tpm[,1])  # = 1

#check colsums after setting -Inf to zero; should total 1e+06
colSums(tmmtpm) #all 1e06
colSums(tpm)  ##all 1e06

#just for fun; how different are tpm and fpkm
plot(log2(fpkm[,1]), log2(tpm[,1]))
#other comparisons
plot(logfpkm_tmm[,1], logtpm[,1])
cor(logfpkm_tmm[,1], logfpkm[,1])
fit <- lm(logfpkm_tmm[,1] ~ logtpm[,1])
summary(fit)
plot(logfpkm_tmm[,1], logfpkm[,1])
cor(logfpkm_tmm[,1], logfpkm[,1])
fit <- lm(logfpkm_tmm[,1] ~ logfpkm[,1])
summary(fit)
#Do the same tests with an xpress count matrix

###################################

#Xpress DGEobj
rm(list=ls())
dge_x <- readRDS("../X20200Gene_dgeobj.RDS")
#filter out genelength=0 genes
el <- getItem(dge_x, "effectiveLength")
el <- rowMeans(el)
idx <- el == 0
print(glue::glue("zerolength gene count = {sum(idx)}."))
dim(dge_x)
dge_zerolength <- dge_x[idx,]
dge_x <- dge_x[!idx,]
dim(dge_x)

#look at counts for zerolength genes
hist(getItem(dge_zerolength, "counts"))
max(getItem(dge_zerolength, "counts"))
summary(getItem(dge_zerolength, "counts")[,1])
View(getItem(dge_zerolength, "counts"))
#So vast majority of zerolenght genes have zero counts

#Now look at the non-zero genelengths
counts <- getItem(dge_x, "counts")
geneData <- getItem(dge_x, "geneData")
design <- getItem(dge_x, "design")
tpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize="tmm", debug=T)
logtpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize="tmm", log=TRUE, debug=T)
fpkm <- convertCounts(counts, unit="fpkm", geneLength=geneData$ExonLength, normalize="none", debug=T)
fpkm_tmm  <- convertCounts(counts, unit="fpkm", geneLength=geneData$ExonLength, normalize="tmm", debug=T)
logfpkm <- convertCounts(counts, unit="fpkm", geneLength=geneData$ExonLength, normalize="none", log=TRUE, debug=T)
logfpkm_tmm <- convertCounts(counts, unit="fpkm", geneLength=geneData$ExonLength, normalize="tmm", log=TRUE, debug=T)

inspect(counts)
inspect(tpm)
inspect(logtpm)
inspect(fpkm)
inspect(fpkm_tmm)
inspect(logfpkm)
inspect(logfpkm_tmm)
#zeros go away for logged data because of the prior.

#the tpm.direct function calculates tpm by the equation without help from any edgeR functions
effectiveLength <- getItem(dge_x, "effectiveLength")
tpm_direct <- tpm.direct(counts, genelength = effectiveLength, collapse=T)
inspect(tpm_direct)
plot(log2(tpm[,1]), log2(tpm_direct[,1]))
diff <- tpm[,1] - tpm_direct[,1]
hist(diff)
max(diff)
#independently calculated TPMs agree to 2e-13

################################

#compare to Xpress TPM data

##get TPM from Xpress
library(rXpress)
projectNumber <- 20200
variables <- rXpress::getVariables(projectNumber)
#which genemodel
variables$projectInfo$DESIGNS
#which attributes
variables$projectInfo$ATTRIBUTES
#241=isoformPct, 203=fpkm, 201=tpm, 205=counts, 207=normalizedCounts, 227=effectiveLength
attributes <- getAttributes(projectNumber)
#retrieve counts (attribute 201 TPM)
xpressData <- rXpress::getXpressData(projectNumber, attribute=201)
varSetNames <- rXpress::getVarSetNames(projectNumber)
Level <- varSetNames[1]  #1 for gene, 2 for isoforms
# MyRowData <- xpressData$varSets[[Level]]$z
# MyColData <- xpressData$varSets[[Level]]$x
TPM_xpress = xpressData$varSets[[Level]]$y

# Remove zero length genes
dge_x <- readRDS("../X20200Gene_dgeobj.RDS")
dim(TPM_xpress)
dim(dge_x)

effectiveLength <- getItem(dge_x, "effectiveLength")
el <- rowMeans(effectiveLength)
idx <- el == 0
sum(idx)

#remove zerolength genes from both sets
TPM_xpress <- TPM_xpress[!idx,]
dge_x_zerolengenes <- dge_x[idx,]
dge_x <- dge_x[!idx,]
dim(TPM_xpress)
dim(dge_x)
#same length? 

#filter TPM_xpress to only gene in common with dge_x)
idx <- rownames(TPM_xpress) %in% rownames(dge_x)
sum(idx)
TPM_xpress <- TPM_xpress[idx,]
dim(TPM_xpress)
dim(dge_x)
#check for same order # confirm sort and compare.
all(rownames(TPM_xpress) == rownames(dge_x))

logtpm <- convertCounts(getItem(dge_x, "counts"), unit="tpm", 
                     geneLength=getItem(dge_x, "geneData")$ExonLength, 
                     normalize="none", log=TRUE, debug=T, prior.count=0.25)
all(rownames(TPM_xpress) == rownames(logtpm))
plot(log2(TPM_xpress[,1]), logtpm[,1])

#compare TPM_xpress to tpm.direct next