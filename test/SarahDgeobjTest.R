dgeobj <- readRDS("X20200Gene_dgeobj.RDS")
library(DGEobj)
library(DGE.Tools2)
library(JRTutil)

counts = getItem(dgeobj, 'counts')
genelength = getItem(dgeobj, 'effectiveLength')
fpk = convertCounts(counts, unit = 'fpk', geneLength = genelength)
logtpm = convertCounts(counts, unit = 'TPM', geneLength = genelength, log = T, normalize = 'TMM')
inspect(counts)
inspect(fpk)
#count and FPK look the same  (note:  ZC is zero count and ZF is zero fraction)
#counts and fpk are fully populated but have but have ~62% zeros.  Zero counts goes to -Inf when you log?
inspect(logtpm)
6766131+877329
#[1] 7643460
#indeed tpm is messed up; mostly -Inf and some NAs
##Let's try not logging
tpm <- convertCounts(counts, unit = 'TPM', geneLength = genelength, log = F, normalize = 'TMM')
inspect(tpm)
#I'd expect tpm (not logged) to look like counts or fpk
#At least one reason why it doesn't is that many of the genelengths are zero
#zero counts divided by anything is 0 and log2(0 = -Inf)
inspect(genelength)
#The rowmeans of genelength is used in the calculation
genelengthVector <- rowMeans(genelength, na.rm = TRUE)
inspect(genelengthVector)

#so Anything with 0 counts or 0 genelength will go to -Inf or Inf  
# zero divided by 0 is NaN
# But cells that are numeric for both genelength and counts should yield a real values)
# #But there is data that is finite for both counts and genelengths that should yield numeric results.
