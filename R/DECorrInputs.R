#' Generate all the input files for TEKRABber downstream analysis
#' @description Generate all the inputs files for differentially expressed 
#' orthologous genes/TEs analysis, and for correlation analysis. The output 
#' is a list containing 6 dataframes. 
#' @usage DECorrInputs(fetchData)
#' @param fetchData output list from TEKRABber::orthologScale()
#' @return create inputs for DE analysis and correlations: 
#' (1) geneInputDESeq2 (2) teInputDESeq2 (3) geneCorrInputRef 
#' (4) geneCorrInputCompare (5) TECorrInputRef (6) TECorrInputCompare
#' @export
#' @importFrom utils write.table
#' @examples
#' data(speciesCounts)
#' data(hg38_panTro6_rmsk)
#' hmGene <- speciesCounts$hmGene
#' chimpGene <- speciesCounts$chimpGene
#' hmTE <- speciesCounts$hmTE
#' chimpTE <- speciesCounts$chimpTE
#' 
#' ## For demonstration, here we only select 1000 rows to save time
#' set.seed(1234)
#' hmGeneSample <- hmGene[sample(nrow(hmGene), 1000), ]
#' chimpGeneSample <- chimpGene[sample(nrow(chimpGene), 1000), ]
#' 
#' fetchData <- orthologScale(
#'     speciesRef = "hsapiens",
#'     speciesCompare = "ptroglodytes",
#'     geneCountRef = hmGeneSample,
#'     geneCountCompare = chimpGeneSample,
#'     teCountRef = hmTE,
#'     teCountCompare = chimpTE,
#'     rmsk = hg38_panTro6_rmsk
#' )
#' 
#' inputBundle <- DECorrInputs(fetchData)
DECorrInputs <- function(fetchData) {
    
    ## save two input for correlation and DE
    i <- c("refEnsemblID", "compareEnsemblID")
    orthologTable_ID <- fetchData$orthologTable[, i]
    
    ## create input for DESeq2
    geneInputDESeq2 <- merge(
        orthologTable_ID, fetchData$geneRef, 
        by = c("refEnsemblID", "refEnsemblID")
    )
    
    geneInputDESeq2 <- merge(
        geneInputDESeq2, fetchData$geneCompare, 
        by = c("compareEnsemblID", "compareEnsemblID")
    )
    
    teInputDESeq2 <- merge(fetchData$teRef, fetchData$teCompare, 
                           by = c("teName", "teName"))
    
    ## remove ensembl ID column
    refCount <- ncol(fetchData$geneRef) - 1
    compareCount <- ncol(fetchData$geneCompare) - 1
    
    ## remove duplicated
    ## rename row names and format DESeq2 input
    geneInputDESeq2 <- 
        geneInputDESeq2[!duplicated(geneInputDESeq2[, "refEnsemblID"]), ]
    rownames(geneInputDESeq2) <- geneInputDESeq2[, "refEnsemblID"]
    geneInputDESeq2 <- geneInputDESeq2[, 3:ncol(geneInputDESeq2)]
    
    teInputDESeq2 <- teInputDESeq2[!duplicated(teInputDESeq2[, "teName"]), ]
    rownames(teInputDESeq2) <- teInputDESeq2[, "teName"]
    teInputDESeq2 <- teInputDESeq2[, 2:ncol(teInputDESeq2)]
    
    output <- list(
        "geneInputDESeq2" = geneInputDESeq2,
        "teInputDESeq2" = teInputDESeq2
    )
    
    output
}
