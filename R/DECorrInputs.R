#' Generate all the input files for TEKRABber downstream analysis
#' @description Generate all the inputs files for differentially expressed 
#' genes/TEs analysis, and for correlation analysis. The output 
#' is a list containing 6 dataframes. 
#' @usage DECorrInputs(orthologTable, scaleFactor, geneCountRef, 
#' geneCountCompare, teCountRef, teCountCompare)
#' @param orthologTable orthologTable output from using orthologScale()
#' @param scaleFactor scaleFactor output from using orthologScale()
#' @param geneCountRef Gene counts from your reference species. First 
#' column should be Ensmebl gene ID.
#' @param geneCountCompare Gene counts from the species you want to compare. 
#' First column should also be Ensembl gene ID.
#' @param teCountRef TE counts from your reference species. First column 
#' should be TE's name.
#' @param teCountCompare TE counts from the species you want to compare. 
#' First column should also be TE's name.
#' @return create inputs for DE analysis and correlations: 
#' (1) geneInputDESeq2 (2) teInputDESeq2 (3) geneCorrInputRef 
#' (4) geneCorrInputCompare (5) TECorrInputRef (6) TECorrInputCompare
#' @export
#' @importFrom utils write.table
#' @examples
#' data(speciesCounts)
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
#'     geneCountCompare = chimpGeneSample
#' )
#' 
#' inputBundle <- DECorrInputs(
#'     orthologTable=fetchData$orthologTable,
#'     scaleFactor=fetchData$scaleFactor,
#'     geneCountRef=hmGene,
#'     geneCountCompare=chimpGene,
#'     teCountRef=hmTE,
#'     teCountCompare=chimpTE
#' )
DECorrInputs <- function(
    orthologTable, scaleFactor, 
    geneCountRef, geneCountCompare, 
    teCountRef, teCountCompare) {
    
    norm_scale <- function(x) {
        return(round(x / scaleFactor))
    }
    
    ## normalize data
    for (i in seq_len(ncol(geneCountCompare))[-1]){
        geneCountCompare[, i] <- norm_scale(geneCountCompare[, i])
    }
    
    for (i in seq_len(ncol(teCountCompare))[-1]){
        teCountCompare[, i] <- norm_scale(teCountCompare[, i])
    }
    
    ## save two input for correlation and DE
    i <- c("refEnsemblID", "compareEnsemblID")
    orthologTable_ID <- orthologTable[, i]
    
    ## create input for DESeq2
    colnames(geneCountRef)[1] <- "refEnsemblID"
    colnames(geneCountCompare)[1] <- "compareEnsemblID"
    colnames(teCountRef)[1] <- "teName"
    colnames(teCountCompare)[1] <- "teName"
    
    geneInputDESeq2 <- merge(
        orthologTable_ID, geneCountRef, 
        by = c("refEnsemblID", "refEnsemblID")
    )
    
    geneInputDESeq2 <- merge(
        geneInputDESeq2, geneCountCompare, 
        by = c("compareEnsemblID", "compareEnsemblID")
    )
    
    teInputDESeq2 <- merge(
        teCountRef, teCountCompare, 
        by = c("teName", "teName")
    )
    
    ## remove ensembl ID column
    refCount <- ncol(geneCountRef) - 1
    compareCount <- ncol(geneCountCompare) - 1
    
    ## remove duplicated
    ## rename row names and format DESeq2 input
    geneInputDESeq2 <- 
        geneInputDESeq2[!duplicated(geneInputDESeq2[, "refEnsemblID"]), ]
    rownames(geneInputDESeq2) <- geneInputDESeq2[, "refEnsemblID"]
    geneInputDESeq2 <- geneInputDESeq2[, 3:ncol(geneInputDESeq2)]
    
    teInputDESeq2 <- teInputDESeq2[!duplicated(teInputDESeq2[, "teName"]), ]
    rownames(teInputDESeq2) <- teInputDESeq2[, "teName"]
    teInputDESeq2 <- teInputDESeq2[, 2:ncol(teInputDESeq2)]
    
    ## in case the number is not integers
    geneInputDESeq2 <- round(geneInputDESeq2)
    teInputDESeq2 <- round(teInputDESeq2)
    
    geneCorrInputRef <- geneInputDESeq2[, seq_len(refCount)]
    geneCorrInputCompare <- 
        geneInputDESeq2[, (refCount + 1):ncol(geneInputDESeq2)]
    
    teCorrInputRef <- teInputDESeq2[, seq_len(refCount)]
    teCorrInputCompare <- 
        teInputDESeq2[, (refCount + 1):ncol(teInputDESeq2)]
    
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