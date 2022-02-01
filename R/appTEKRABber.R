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
#' ## hmchimpDE is generated from using DEgeneTE()
#' ## hmCorrResult and chimpCorrResult are generated from using corrOrthologTE()
#' ## meta is the same metadata you used for DE analysis
#' \donttest{
#' appTEKRABber(
#'   DEresult = hmchimpDE,
#'   corrRef = hmCorrResult,
#'   corrCompare = chimpCorrResult,
#'   metadata = meta
#' )}
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