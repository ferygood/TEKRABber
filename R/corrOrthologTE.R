#' Estimate correlation comparing orthologs and TEs
#' @description To estimate correlation comparing each ortholog to each TE from inputs.
#' You can sepcify the correlation and adjusted p-value methods (see details in parameters).
#' It will also generate a output csv file in ./results directory.
#' a output csv file with the correlation result.
#' @usage corrOrthologTE(geneInput, teInput, corrMethod = "pearson", padjMethod = "fdr", filename)
#' @param geneInput gene count input for correlation from using DECorrInputs()
#' @param teInput te count input for correlation from using DECorrInputs()
#' @param corrMethod correlation method, inlcuding pearson, kendall, spearman. Default is pearson
#' @param padjMethod method to return p-values adjusted, and default is fdr. See ?p.adjust
#' @param filename specify a csv file name such as "correlationResult.csv"
#' @return a dataframe includes Pearson's correlation coefficient, pvalue, padj
#' @export
#' @examples
#' library(SummarizedExperiment)
#' data(speciesCorr)
#' hmGeneCorrInput <- assay_tekcorrset(speciesCorr, "gene", "human")
#' hmTECorrInput <- assay_tekcorrset(speciesCorr, "te", "human")
#'
#' #corrOrthologTE(
#' #    geneInput=hmGeneCorrInput,
#' #    teInput=hmTECorrInput,
#' #    corrMethod="pearson",
#' #    padjMethod="fdr",
#' #    filename="correlationResult.csv"
#' #)
corrOrthologTE <- function(geneInput, teInput, corrMethod = "pearson", padjMethod="fdr", filename){
    dir.create("./results")
    df.ortholog <- t(geneInput)
    df.te <- t(teInput)
    
    df.corr <- rcpp_corr(df.ortholog, df.te, corrMethod)
    colnames(df.corr)[1:2] <- c("geneName", "teName")
    df.corr <- df.corr %>%
        mutate(padj = p.adjust(pvalue, method=padjMethod))
    rownames(df.corr) <- c(1:nrow(df.corr))
    write.table(df.corr, file = paste0("results/", filename), sep=",")
    
    df.corr
    
}