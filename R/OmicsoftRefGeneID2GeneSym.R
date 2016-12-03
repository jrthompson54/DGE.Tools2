#JRT 14Nov2015
#Convert Omicsoft Refgene GeneID to GeneSymbol
### Function OmicsoftRefGeneID2GeneSym ###
#' Function OmicsoftRefGeneID2GeneSym:  Convert Omicsoft RefGene GeneIDs to GeheSym
#'
#' Maps Omicsoft RegGene IDs (e.g. NOC2L_1_-_1) to official gene symbol (e.g. NOC2L).
#' Omicsoft appends chr information to the genesymbol to make them unique.
#' This procedure simply right trims everything past the underscore
#'
#' @author John Thompson, \email{john.thompson@@bms.com}
#' @keywords gene symbol, Entrez, GeneID
#'
#' @param GeneID a vector of Omicsoft GeneIDs
#'
#' @return A vector with GeneSymbols
#'
#' @examples
#' #Senario GeneID is the rowname.  Add a GeneSym column
#' MyDF$MyGeneSym = OmicsoftRefGeneID2GeneSym (rownames(MyDF))
#'
#' @export
OmicsoftRefGeneID2GeneSym <- function(GeneID) {
  gs = strsplit(GeneID, "_" )  #official genesym is first element of each list
  GeneSym = lapply(gs, '[[', 1) %>% unlist
}
