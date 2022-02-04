#' Estimate correlation comparing orthologs and TEs
#' @description To estimate correlation comparing orthologs and TEs one-by-one from inputs.
#' You can specify the correlation and adjusted p-value methods (see details in parameters).
#' If you want to save your outputs instead just return them, please specify the fileDir and 
#' fileName with the extension .csv. The default fileName is TEKRABber_geneTECorrReusult.csv.
#' a output csv file with the correlation result.
#' @usage corrOrthologTE(geneInput, teInput, corrMethod = "pearson", padjMethod = "fdr", fileDir=NULL, fileName="TEKRABber_geneTECorrResult.csv")
#' @param geneInput gene count input for correlation from using DECorrInputs()
#' @param teInput te count input for correlation from using DECorrInputs()
#' @param corrMethod correlation method, including pearson, kendall, spearman. Default is pearson
#' @param padjMethod method to return p-values adjusted, and default is fdr. See ?p.adjust
#' @param fileDir the name and path of directory for saving output files. Default is NULL.
#' @param fileName the name for saving output files. Default is "TEKRABber_geneTECorrResult.csv"
#' @return a dataframe includes correlation coefficient, pvalue, padj
#' @useDynLib TEKRABber
#' @importFrom stats p.adjust
#' @importFrom Rcpp sourceCpp
#' @export
#' 
#' @examples
#' library(SummarizedExperiment)
#' data(speciesCorr)
#' hmGeneCorrInput <- assay_tekcorrset(speciesCorr, "gene", "human")
#' hmTECorrInput <- assay_tekcorrset(speciesCorr, "te", "human")
#'
#' \donttest{
#' corrOrthologTE(
#'     geneInput=hmGeneCorrInput,
#'     teInput=hmTECorrInput,
#'     corrMethod="pearson",
#'     padjMethod="fdr",
#'     fileDir=".",
#'     fileName="correlationResult.csv"
#' )}
corrOrthologTE <- function(
    geneInput, 
    teInput, 
    corrMethod = "pearson", 
    padjMethod="fdr", 
    fileDir=NULL, 
    fileName="TEKRABber_geneTECorrResult.csv"){

    df.ortholog <- t(geneInput)
    df.te <- t(teInput)
    
    df.corr <- rcpp_corr(df.ortholog, df.te, corrMethod)
    colnames(df.corr)[c(1,2)] <- c("geneName", "teName")
    df.corr$padj <- p.adjust(
        df.corr$pvalue,
        method=padjMethod
    )
    rownames(df.corr) <- seq_len(nrow(df.corr))
    
    # if user has specify file directory and file name
    if (!is.null(fileDir)){
        dir.create(fileDir)
        write.table(df.corr, file = file.path(fileDir, fileName), sep=",")
    }
    
    df.corr
    
}