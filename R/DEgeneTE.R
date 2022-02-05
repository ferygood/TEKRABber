#' Estimate differentially expressed genes and TEs
#' @description To estimate differentially expressed genes and TEs, DEgeneTE() takes
#' gene inputs and TE inputs from the results using the DECorrInputs function. You need to
#' specify your metadata, contrastVector, and expDesign based on your design. If you 
#' also want to save the output, please specify the fileDir parameter.
#' @usage DEgeneTE(geneTable, teTable, metadata, contrastVector, expDesign=TRUE, fileDir=NULL)
#' @param geneTable gene input table from using DECorrInputs()
#' @param teTable TE input table from using DECorrInputs()
#' @param metadata a one column dataframe with rownames same as the column name of gene/te count table. Column name must be \strong{species} or \strong{experiment}.
#' @param contrastVector your experiment design, i.e. c("species", "human", "chimpanzee")
#' @param expDesign Logic value for comparing between or within species. \strong{TRUE} for comparing between two species, and \strong{FALSE} for comparing between control and treatment.
#' @param fileDir the name and path of directory for saving output files. Default is NULL.
#' @return return DESeq2 res and normalized gene counts.
#' @import apeglm
#' @export
#' @examples
#' ## comparing between species: 
#' ## (1) set expDesign = TRUE (2) column name of metadata needs to be "species".
#' 
#' data(speciesCounts)
#' hmGene <- speciesCounts$hmGene
#' hmTE <- speciesCounts$hmTE
#' chimpGene <- speciesCounts$chimpGene
#' chimpTE <- speciesCounts$chimpTE
#' 
#' data(fetchDataHmChimp)
#' fetchData <- fetchDataHmChimp
#' 
#' inputBundle <- DECorrInputs(
#'   orthologTable = fetchData$orthologTable,
#'   scaleFactor = fetchData$scaleFactor,
#'   geneCountRef = hmGene,
#'   geneCountCompare = chimpGene,
#'   teCountRef = hmTE,
#'   teCountCompare = chimpTE
#' )
#' 
#' meta <- data.frame(species=c(rep("human", ncol(hmGene) - 1), 
#'                              rep("chimpanzee", ncol(chimpGene) - 1))
#' )
#' rownames(meta) <- colnames(inputBundle$geneInputDESeq2)
#' meta$species <- factor(meta$species, levels = c("human", "chimpanzee"))
#' 
#' hmchimpDE <- DEgeneTE(
#'   geneTable = inputBundle$geneInputDESeq2,
#'   teTable = inputBundle$teInputDESeq2,
#'   metadata = meta,
#'   contrastVector = c("species", "human", "chimpanzee"),
#'   expDesign = TRUE
#' )
DEgeneTE <- function(geneTable, teTable, metadata, contrastVector, expDesign = TRUE, fileDir = NULL) {
    deseq2 <- function(cts, coldata) {
        if (all(rownames(coldata) == colnames(cts))) {
            # round counts to integers
            cts <- round(cts) 
            
            dds <- c()
            if (expDesign == TRUE) {
               dds <- DESeq2::DESeqDataSetFromMatrix(
                   countData = cts, 
                   colData = coldata, 
                   design = ~species
               ) 
            } else if (expDesign == FALSE) {
                dds <- DESeq2::DESeqDataSetFromMatrix(
                    countData = cts,
                    colData = coldata,
                    design = ~experiment
                )
            }
            
            ## pre-filter and normalized
            keep <- rowSums(DESeq2::counts(dds)) >= 10
            dds <- dds[keep, ]
            dds_dataset <- dds
            dds <- DESeq2::DESeq(dds)
            normalized_counts <- DESeq2::counts(dds, normalized = TRUE)
            
            coefName <- DESeq2::resultsNames(dds)[2]
            res <- DESeq2::results(dds, contrast = contrastVector, alpha = 0.05)
            res <- DESeq2::lfcShrink(dds, coef = coefName, type = "apeglm")
            
            result <- list("dds" = dds_dataset, "normalized_counts" = normalized_counts, "res" = res)
            result
        }
    }
    
    ## run analysis
    
    geneDE <- deseq2(geneTable, metadata)
    teDE <- deseq2(teTable, metadata)
    
    ## save files if directory is specified
    if (!is.null(fileDir)){
        dir.create(fileDir)
        write.table(data.frame(geneDE$normalized_counts), file = file.path(fileDir, "geneDESeq2Log2.csv", sep = ","))
        write.table(data.frame(geneDE$res), file = file.path(fileDir, "geneDESeq2results.csv", sep = ","))
        write.table(data.frame(teDE$normalized_counts), file = file.path(fileDir, "teDESeq2Log2.csv", sep = ","))
        write.table(data.frame(teDE$res), file = file.path(fileDir, "teDESeq2results.csv", sep = ","))
    }
    
    output <- list(
        "gene_dds" = geneDE$dds,
        "gene_res" = geneDE$res,
        "normalized_gene_counts" = geneDE$normalized_counts,
        "te_dds" = teDE$dds,
        "te_res" = teDE$res,
        "normalized_te_counts" = teDE$normalized_counts
    )
    
    output
    
}