#' Generate all the inputs for TEKRABber downstream analysis
#' @description Generate inputs files for differentially expressed genes/TEs analysis,
#' and for correlation analysis. The output is a list contain all input files. If you
#' prefer using gene name instead of Ensembl gene ID, you can convert them before running
#' this function.
#' @usage DECorrInputs(orthologTable, scaleFactor, geneCountRef, geneCountCompare, teCountRef, teCountCompare)
#' @param orthologTable orthologTable output from using orthologScale()
#' @param scaleFactor scaleFactor output from using orthologScale()
#' @param geneCountRef Gene counts from your reference species. First column should be Ensmebl gene ID
#' @param geneCountCompare Gene counts from the species you want to compare. First column should also be Ensembl gene ID
#' @param teCountRef TE counts from your reference species. First column should be TE's name
#' @param teCountCompare TE counts from the species you want to compare. First column should also be TE's name
#' @return create inputs for DE analysis and correlations: (1) geneInputDESeq2 (2) teInputDESeq2 (3) geneCorrInputRef (4) geneCorrInputCompare (5) TECorrInputRef (6) TECorrInputCompare
#' @export
#' @importFrom magrittr %>%
#' @importFrom dplyr mutate inner_join across
#' @importFrom utils write.table
#' @examples
#' data(speciesCounts)
#' 
#' hmGene <- speciesCounts$hmGene
#' chimpGene <- speciesCounts$chimpGene
#' chimpGene <- speciesCounts$chimpGene
#' chimpTE <- speciesCounts$chimpTE
#' \donttest{
#' fetchData <- orthologScale(
#'   speciesRef = "hsapiens",
#'   speciesCompare = "ptroglodytes",
#'   geneCountRef = hmGene,
#'   geneCountCompare = chimpGene
#' )}
#' \donttest{
#' inputBundle <- DECorrInputs(
#'   orthologTable=fetchData$orthologTable,
#'   scaleFactor=fetchData$scaleFactor,
#'   geneCountRef=hmGene,
#'   geneCountCompare=chimpGene,
#'   teCountRef=hmTE,
#'   teCountCompare=chimpTE
#' )}
DECorrInputs <- function(orthologTable, scaleFactor, geneCountRef, 
                         geneCountCompare, teCountRef, teCountCompare) {
    norm_scale <- function(x) {
        return(round(x / scaleFactor))
    }
    
    geneCountCompare <- geneCountCompare %>%
        mutate(across(2:ncol(geneCountCompare), norm_scale))
    
    teCountCompare <- teCountCompare %>%
        mutate(across(2:ncol(teCountCompare), norm_scale))
    
    ## save two input for correlation and DE
    orthologTable_ID <- orthologTable[, c(3,7)]
    
    ## create input for DESeq2
    colnames(geneCountRef)[1] <- "refEnsemblID"
    colnames(geneCountCompare)[1] <- "compareEnsemblID"
    colnames(teCountRef)[1] <- "teName"
    colnames(teCountCompare)[1] <- "teName"
    
    geneInputDESeq2 <- inner_join(orthologTable_ID, geneCountRef, by = "refEnsemblID")
    geneInputDESeq2 <- inner_join(geneInputDESeq2, geneCountCompare, by = "compareEnsemblID")
    
    teInputDESeq2 <- inner_join(teCountRef, teCountCompare, by = "teName")
    
    ## remove duplicated rows in case
    refCount <- ncol(geneCountRef[, 2:ncol(geneCountRef)])
    compareCount <- ncol(geneCountCompare[, 2:ncol(geneCountCompare)])
    
    ## rename row names and format DESeq2 input
    geneInputDESeq2 <- geneInputDESeq2[!duplicated(geneInputDESeq2[, 1]), ]
    rownames(geneInputDESeq2) <- geneInputDESeq2[, 1]
    geneInputDESeq2 <- geneInputDESeq2[, 3:ncol(geneInputDESeq2)]
    
    teInputDESeq2 <- teInputDESeq2[!duplicated(teInputDESeq2[, 1]), ]
    rownames(teInputDESeq2) <- teInputDESeq2[, 1]
    teInputDESeq2 <- teInputDESeq2[, 2:ncol(teInputDESeq2)]
    
    ## in case the number is not integers
    geneInputDESeq2 <- round(geneInputDESeq2)
    teInputDESeq2 <- round(teInputDESeq2)
    
    #write.table(geneInputDESeq2, file = "results/geneInputDESeq2.csv", sep = ",")
    #write.table(teInputDESeq2, file = "results/teInputDESeq2.csv", sep = ",")
    
    geneCorrInputRef <- geneInputDESeq2[, 1:refCount]
    geneCorrInputCompare <- geneInputDESeq2[, (refCount + 1):ncol(geneInputDESeq2)]
    teCorrInputRef <- teInputDESeq2[, 1:refCount]
    teCorrInputCompare <- teInputDESeq2[, (refCount + 1):ncol(teInputDESeq2)]
    
    output <- list(
        "geneInputDESeq2" = geneInputDESeq2,
        "teInputDESeq2" = teInputDESeq2,
        "geneCorrInputRef" = geneCorrInputRef,
        "geneCorrInputCompare" = geneCorrInputCompare,
        "TECorrInputRef" = teCorrInputRef,
        "TECorrInputCompare" = teCorrInputCompare
    )
    
    output
}