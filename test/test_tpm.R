library(DGEobj)
library(DGE.Tools2)
library(JRTutil)

#omicsoft DGEobj
dge_as <- readRDS("../DGEobj.RDS")
counts <- getItem(dge_as, "counts")
geneData <- getItem(dge_as, "geneData")
design <- getItem(dge_as, "design")
tpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize=TRUE)
logtpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize=TRUE, log=TRUE)
fpkm <- convertCounts(counts, unit="fpkm", geneLength=geneData$ExonLength, normalize=TRUE)
inspect(counts)
inspect(tpm)
inspect(logtpm)
hist(logtpm[,1], breaks=100)

plot(logtpm[,1], log2(tpm[,1]))
# > inspect(counts)
# Item    Value
# 1           Class   matrix
# 2          Length   260964
# 3             Dim 14498:18
# 4  complete.cases    14498
# 5              ZC     1441
# 6              ZF    0.006
# 7             Inf        0
# 8          NegInf        0
# 9              NA        0
# 10            NAN        0
# > inspect(tpm)
# Item    Value
# 1           Class   matrix
# 2          Length   260964
# 3             Dim 14498:18
# 4  complete.cases    14498
# 5              ZC     1441
# 6              ZF    0.006
# 7             Inf        0
# 8          NegInf        0
# 9              NA        0
# 10            NAN        0
# > inspect(logtpm)
# Item    Value
# 1           Class   matrix
# 2          Length   260964
# 3             Dim 14498:18
# 4  complete.cases    14498
# 5              ZC        0
# 6              ZF        0
# 7             Inf        0
# 8          NegInf     1441
# 9              NA        0
# 10            NAN        0
# Note that Zero counts converts  -Inf logtpm  as it should

#Confirm that TPM is the same as TMM TPM
tmmtpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize=TRUE)
tpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize="none")
plot(log2(tmmtpm[,1]), log2(tpm[,1]))
diff <- tmmtpm - tpm
hist(diff[,1]) #not all zeros but very close.

#check colsums after setting -Inf to zero; should total 1e+06
idx <- tpm == -Inf
tpm[idx] <- 0
colSums(tpm)

idx <- tmmtpm == -Inf
tmmtpm[idx] <- 0
colSums(tmmtpm)

#Do the same tests with an xpress count matrix
#
#Xpress DGEobj
dge_x <- readRDS("../X20200Gene_dgeobj.RDS")
counts <- getItem(dge_x, "counts")
geneData <- getItem(dge_x, "geneData")
design <- getItem(dge_x, "design")
effectiveLength <- getItem(dge_x, "effectiveLength")
tpm1 <- convertCounts(counts, unit="tpm", geneLength=effectiveLength, normalize=TRUE)
tpm2 <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize=TRUE)
logtpm <- convertCounts(counts, unit="tpm", geneLength=geneData$ExonLength, normalize=TRUE, log=TRUE)
inspect(counts)
inspect(tpm1)
inspect(logtpm)

tpmc <- tpm.classic(counts, genelength = effectiveLength, collapse=F)
inspect(tpmc)
