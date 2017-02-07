#test convertCounts

counts <- assay(RSE, "Counts")
genelength <- mcols(RSE)$ExonLength

source('~/R/lib/pkgsrc/DGE.Tools2/R/convertCounts.R')

fkpm <- convertCounts(counts, "FPKM", geneLength=genelength)
logfpkm <- convertCounts(counts, "FPKM", geneLength=genelength, log=TRUE)
tpm <- convertCounts(counts, "TPM", geneLength=genelength)
logtpm <- convertCounts(counts, "TPM", geneLength=genelength, log=TRUE)
cpm <- convertCounts(counts, "CPM", geneLength=genelength)
logcpm <- convertCounts(counts, "CPM", geneLength=genelength, log=TRUE)
fpk <- convertCounts(counts, "FPK", geneLength=genelength)
logfpk <- convertCounts(counts, "FPK", geneLength=genelength, log=TRUE)
zfkpm <- convertCounts(counts, "ZFPKM", geneLength=genelength)
logzfpkm <- convertCounts(counts, "ZFPKM", geneLength=genelength, log=TRUE)

fkpm <- convertCounts(counts, "FPKM", geneLength=genelength, normalize=TRUE)
logfpkm <- convertCounts(counts, "FPKM", geneLength=genelength, log=TRUE, normalize=TRUE)
tpm <- convertCounts(counts, "TPM", geneLength=genelength, normalize=TRUE)
logtpm <- convertCounts(counts, "TPM", geneLength=genelength, log=TRUE, normalize=TRUE)
cpm <- convertCounts(counts, "CPM", geneLength=genelength, normalize=TRUE)
logcpm <- convertCounts(counts, "CPM", geneLength=genelength, log=TRUE, normalize=TRUE)
fpk <- convertCounts(counts, "FPK", geneLength=genelength, normalize=TRUE)
logfpk <- convertCounts(counts, "FPK", geneLength=genelength, log=TRUE, normalize=TRUE)
zfkpm <- convertCounts(counts, "ZFPKM", geneLength=genelength, normalize=TRUE)
logzfpkm <- convertCounts(counts, "ZFPKM", geneLength=genelength, log=TRUE, normalize=TRUE)


fkpm <- convertCounts(counts, "FPKM", geneLength=genelength, normalize="RLE")
logfpkm <- convertCounts(counts, "FPKM", geneLength=genelength, log=TRUE, normalize="RLE")
tpm <- convertCounts(counts, "TPM", geneLength=genelength, normalize="RLE")
logtpm <- convertCounts(counts, "TPM", geneLength=genelength, log=T, normalize="RLE")
cpm <- convertCounts(counts, "CPM", geneLength=genelength, normalize="RLE")
logcpm <- convertCounts(counts, "CPM", geneLength=genelength, log=T, normalize="RLE")
fpk <- convertCounts(counts, "FPK", geneLength=genelength, normalize="RLE")
logfpk <- convertCounts(counts, "FPK", geneLength=genelength, log=TRUE, normalize="RLE")
zfkpm <- convertCounts(counts, "ZFPKM", geneLength=genelength, normalize="RLE")
logzfpkm <- convertCounts(counts, "ZFPKM", geneLength=genelength, log=TRUE, normalize="RLE")
