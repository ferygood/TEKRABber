#' Visualize TEKRABber results with shiny app
#' @description To help user explore their results using TEKRABber, it visualize the results
#' using a self-written shiny app with two tabs, including the expression and correlation of
#' genes and TEs.
#' @usage appTEKRABber(DEresult, corrRef, corrCompare, metadata)
#' @param DEresult the output variable from using DEgeneTE()
#' @param corrRef the correlation result of your reference species using corrOthologTE()
#' @param corrCompare the correlation result of your compare species using corrOrthologTE()
#' @param metadata the same metadata you use for DEgeneTE()
#' 
#' @importFrom shiny runApp
#' @export
#' 
#' @return an app can display differentially expressed genes/TE and the correlation results
#' @examples 
#' data(ctInputDE)
#' geneInputDE <- ctInputDE$gene
#' teInputDE <- ctInputDE$te
#' 
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
#' 
#' #library(SummarizedExperiment)
#' #data(ctCorr)
#' #geneConCorrInput <- assay_tekcorrset(ctCorr, "gene", "control")
#' #teConCorrInput <- assay_tekcorrset(ctCorr, "te", "control")
#' #geneTreatCorrInput <- assay_tekcorrset(ctCorr, "gene", "treatment")
#' #teTreatCorrInput <- assay_tekcorrset(ctCorr, "te", "treatment")
#' 
#' #controlCorr <- corrOrthologTE(
#' #  geneInput = geneConCorrInput,
#' #  teInput = teConCorrInput,
#' #  corrMethod = "pearson",
#' #  padjMethod = "fdr",
#' #  filename = "controlCorrResult.csv"
#' #)
#' 
#' #treatmentCorr <- corrOrthologTE(
#' #  geneInput = geneTreatCorrInput,
#' #  teInput = teTreatCorrInput,
#' #  corrMethod = "pearson",
#' #  padjMethod = "fdr",
#' #  filename = "treatmentCorrResult.csv"
#' #)
#' 
#' #appTEKRABber(
#' #  DEresult = resultDE,
#' #  corrRef = controlCorr,
#' #  corrCompare = treatmentCorr,
#' #  metadata = metaExp
#' #)
#' 
appTEKRABber <- function(DEresult, corrRef, corrCompare, metadata) {
    geneInput <- DEresult$normalized_gene_counts
    geneRes <- DEresult$gene_res
    teInput <- DEresult$normalized_te_counts
    teRes <- DEresult$te_res
    
    rownames(geneInput) <- gsub("[()]", "", rownames(geneInput))
    rownames(teInput) <- gsub("[()]", "", rownames(teInput))
    rownames(geneRes) <- gsub("[()]", "", rownames(geneRes))
    rownames(teRes) <- gsub("[()]", "", rownames(teRes))
    corrRef$geneName <- gsub("[()]", "", corrRef$geneName)
    corrRef$teName <- gsub("[()]", "", corrRef$teName)
    corrCompare$geneName <- gsub("[()]", "", corrCompare$geneName)
    corrCompare$teName <- gsub("[()]", "", corrCompare$teName)
    
    # for expression tab: avoid subscript out of bounds with empty expression
    geneChoices <- intersect(corrRef$geneName, rownames(geneInput))
    geneChoices <- intersect(corrCompare$geneName, geneChoices)
    teChoices <- intersect(corrRef$teName, rownames(teInput))
    teChoices <- intersect(corrCompare$teName, teChoices)
    
    runApp(
        appDir=system.file("shinyGUI", package="TEKRABber"),
        launch.browser=TRUE
    )
    
    
}