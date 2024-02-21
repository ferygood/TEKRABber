#' Estimate correlation comparing orthologs and TEs
#' @description To estimate correlation comparing orthologs and TEs one-by-one 
#' from inputs. You can specify the correlation and adjusted p-value methods 
#' (see details in parameters). If you want to save your outputs instead of 
#' just returning them, please specify the fileDir and fileName with the 
#' extension .csv. The default fileName is TEKRABber_geneTECorrReusult.csv.
#' @usage corrOrthologTE(geneInput, teInput, corrMethod = "pearson", 
#' padjMethod = "fdr", fileDir=NULL, fileName="TEKRABber_geneTECorrResult.csv")
#' @param geneInput gene count input for correlation from using DECorrInputs()
#' @param teInput te count input for correlation from using DECorrInputs()
#' @param corrMethod correlation method, including pearson, kendall, spearman. 
#' Default is pearson.
#' @param padjMethod method to return adjusted p-value, and default is fdr. 
#' See ?p.adjust
#' @param numCore number of cores to run parallel. Default is 1. You can use 
#' detectCores() to get how many cores you can use.
#' @param fileDir the name of directory for saving output files. 
#' Default is NULL.
#' @param fileName the name for saving output files. 
#' Default is "TEKRABber_geneTECorrResult.csv"
#' @return a dataframe includes correlation coefficient, pvalue, padj
#' @useDynLib TEKRABber
#' @importFrom stats p.adjust
#' @importFrom Rcpp sourceCpp
#' @importFrom doParallel registerDoParallel stopImplicitCluster
#' @import foreach
#' @export
#' 
#' @examples
#' data(ctInputDE)
#' geneInputDE <- ctInputDE$gene
#' teInputDE <- ctInputDE$te
#' 
#' metaExp <- data.frame(experiment = c(rep("control", 5), rep("treatment", 5)))
#' rownames(metaExp) <- colnames(geneInputDE)
#' metaExp$experiment <- factor(
#'     metaExp$experiment, 
#'     levels = c("control", "treatment")
#' )
#' 
#' resultDE <- DEgeneTE(
#'     geneTable = geneInputDE,
#'     teTable = teInputDE,
#'     metadata = metaExp,
#'     expDesign = FALSE
#' )
#' 
#' controlCorr <- corrOrthologTE(
#'     geneInput = resultDE$geneCorrInputRef[c(1:10),],
#'     teInput = resultDE$teCorrInputRef[c(1:10),],
#'     corrMethod = "pearson",
#'     padjMethod = "fdr"
#' )
#' 
corrOrthologTE <- function(
    geneInput, 
    teInput, 
    corrMethod = "pearson", 
    padjMethod="fdr",
    numCore=1,
    fileDir=NULL, 
    fileName="TEKRABber_geneTECorrResult.csv"){

    # register backend
    registerDoParallel(numCore)
    
    # create empty dataframe
    df.corr <- data.frame(matrix(ncol=4, nrow=0))
    colnames(result_df) <- c("geneName", "teName", "coef", "pvalue")
    
    # calculate correlation with parallel
    df.corr <- foreach(i=1:nrow(geneInput), .combine=rbind) %dopar% {
        foreach(j=1:nrow(teInput), .combine=rbind) %dopar% {
            num1 <- as.numeric(geneInput[i,])
            num2 <- as.numeric(teInput[j,])
            corr <- cor.test(num1, num2, method=corrMethod)
            rbind(
                data.frame(geneName=rownames(geneInput)[i], 
                           teName=rownames(teInput)[j], 
                           coef=corr$estimate, 
                           pvalue=corr$p.value),
                df.corr
            )
        }
    }
    
    stopImplicitCluster()
    
    # calculate p-adjusted value
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