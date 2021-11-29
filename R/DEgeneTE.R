#' Estimate differentially expressed genes and TEs
#' @description To estimate differentially expressed genes and TEs, DEgeneTE() takes
#' gene input and TE input from the result using DECorrInputs function and return results
#' in ./results directory which is automatically generate by this function. You need to
#' specify metadata, contrastVector, and expDesign. See details in example.
#' @usage DEgeneTE(geneTable, teTable, metadata, contrastVector, expDesign = TRUE)
#' @param geneTable gene input table from using DECorrInputs()
#' @param teTable TE input table from using DECorrInputs()
#' @param metadata a one column dataframe with rownames same as the column name of gene/te count table. Column name must be species or experiment.
#' @param contrastVector your experiment design, i.e. c("species", "human", "chimpanzee")
#' @param expDesign Logic value for comparing between or within species. TRUE for comparing between two species, and False for comparing between control and treatment.
#' @return output DESeq2 res and normalized gene counts result in ./results directory
#' @export
#' @examples
#' # if you are comparing between species:
#' # (1) you need to specify the column name to species in the metadata
#' # (2) you can use the output from DECorrInputs()
#' # (3) you need to set "expDesign = TRUE" (see vignettes/TEKRABber.Rmd for details)
#' data(geneInputDE)
#' data(teInputDE)
#' metaExp <- data.frame(experiment = c(rep("control", 5), rep("treatment", 5)))
#' rownames(metaExp) <- colnames(geneInputDE)
#' metaExp$experiment <- factor(metaExp$experiment, levels = c("control", "treatment"))
#'
#' resultDE <- DEgeneTE(
#'   geneTable = geneInputDE,
#'   teTable = teInputDE,
#'   metadata = metaExp,
#'   contrastVector = c("experiment", "control", "treatment"),
#'   expDesign = FALSE
#' )
DEgeneTE <- function(geneTable, teTable, metadata, contrastVector, expDesign = TRUE) {
    deseq2 <- function(cts, coldata) {
        if (all(rownames(coldata) == colnames(cts))) {
            dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~species)
            
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
    
    deseq2_one <- function(cts, coldata) {
        if (all(rownames(coldata) == colnames(cts))) {
            ## round counts to integers
            cts <- round(cts)
            dds <- DESeq2::DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~experiment)
            
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
    if (expDesign == TRUE) {
        geneDE <- deseq2(geneTable, metadata)
        teDE <- deseq2(teTable, metadata)
        
        ## save files
        write.table(data.frame(geneDE$normalized_counts), file = "results/geneDESeq2Log2.csv", sep = ",")
        write.table(data.frame(geneDE$res), file = "results/geneDESeq2results.csv", sep = ",")
        write.table(data.frame(teDE$normalized_counts), file = "results/teDESeq2Log2.csv", sep = ",")
        write.table(data.frame(teDE$res), file = "results/teDESeq2results.csv", sep = ",")
        
        output <- list(
            "gene_dds" = geneDE$dds,
            "gene_res" = geneDE$res,
            "normalized_gene_counts" = geneDE$normalized_counts,
            "te_dds" = teDE$dds,
            "te_res" = teDE$res,
            "normalized_te_counts" = teDE$normalized_counts
        )
        return(output)
    } else if (expDesign == FALSE) {
        ## create results directory if user use this function directly
        dir.create("./results")
        geneDE <- deseq2_one(geneTable, metadata)
        teDE <- deseq2_one(teTable, metadata)
        
        ## save files
        write.table(data.frame(geneDE$normalized_counts), file = "results/geneDESeq2Log2.csv", sep = ",")
        write.table(data.frame(geneDE$res), file = "results/geneDESeq2results.csv", sep = ",")
        write.table(data.frame(teDE$normalized_counts), file = "results/teDESeq2Log2.csv", sep = ",")
        write.table(data.frame(teDE$res), file = "results/teDESeq2results.csv", sep = ",")
        
        output <- list(
            "gene_dds" = geneDE$dds,
            "gene_res" = geneDE$res,
            "normalized_gene_counts" = geneDE$normalized_counts,
            "te_dds" = teDE$dds,
            "te_res" = teDE$res,
            "normalized_te_counts" = teDE$normalized_counts
        )
        return(output)
    } else {
        (
            stop("expDesign should be specified TRUE for comparing speceis, or FALSE for comparing control and experiment.")
        )
    }
}