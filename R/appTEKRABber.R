#' Visualize TEKRABber results with shiny app
#' @description To help user explore their results using TEKRABber, it 
#' visualizes the results using a self-written shiny app with two tabs, 
#' including the expression and correlation of genes and TEs. This function 
#' will create global app-prefix variables to run the app.
#' @usage appTEKRABber(DEresult, corrRef, corrCompare, metadata)
#' @param DEresult the output variable from using DEgeneTE()
#' @param corrRef the correlation result of your reference species using 
#' corrOthologTE()
#' @param corrCompare the correlation result of your compare species using 
#' corrOrthologTE()
#' @param metadata the same metadata you use for DEgeneTE()
#' 
#' @export
#' @return An app to display differentially expressed genes/TEs and the 
#' correlation results
#' @examples
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
#'     orthologTable = fetchData$orthologTable,
#'     scaleFactor = fetchData$scaleFactor,
#'     geneCountRef = hmGene,
#'     geneCountCompare = chimpGene,
#'     teCountRef = hmTE,
#'     teCountCompare = chimpTE
#' )
#' 
#' meta <- data.frame(species=c(rep("human", ncol(hmGene) - 1), 
#'     rep("chimpanzee", ncol(chimpGene) - 1))
#' )
#' rownames(meta) <- colnames(inputBundle$geneInputDESeq2)
#' meta$species <- factor(meta$species, levels = c("human", "chimpanzee"))
#' 
#' hmchimpDE <- DEgeneTE(
#'     geneTable = inputBundle$geneInputDESeq2,
#'     teTable = inputBundle$teInputDESeq2,
#'     metadata = meta,
#'     contrastVector = c("species", "human", "chimpanzee"),
#'     expDesign = TRUE
#' )
#' 
#' data(speciesCorr)
#' hmGeneCorrInput <- assay_tekcorrset(speciesCorr, "gene", "human")
#' hmTECorrInput <- assay_tekcorrset(speciesCorr, "te", "human")
#' chimpGeneCorrInput <- assay_tekcorrset(speciesCorr, "gene", "chimpanzee")
#' chimpTECorrInput <- assay_tekcorrset(speciesCorr, "te", "chimpanzee")

#' hmCorrResult <- corrOrthologTE(
#'     geneInput = hmGeneCorrInput,
#'     teInput = hmTECorrInput,
#'     corrMethod = "pearson",
#'     padjMethod = "fdr"
#' )
# 
#' chimpCorrResult <- corrOrthologTE(
#'     geneInput = chimpGeneCorrInput,
#'     teInput = chimpTECorrInput,
#'     corrMethod = "pearson",
#'     padjMethod = "fdr"
#' )
#' 
#' if (interactive()){
#'     appTEKRABber(
#'         DEresult = hmchimpDE,
#'         corrRef = hmCorrResult,
#'         corrCompare = chimpCorrResult,
#'         metadata = meta
#'     )
#' }
#' 
appTEKRABber <- function(DEresult, corrRef, corrCompare, metadata) {
    
    # create global variables for app-use
    assign("appDE", DEresult, envir = .GlobalEnv)
    assign("appRef", corrRef, envir = .GlobalEnv)
    assign("appCompare", corrCompare, envir = .GlobalEnv)
    assign("appMeta", metadata, envir = .GlobalEnv)
    
    # run shiny app
    shiny::runApp(
        appDir = system.file("shinyGUI", package="TEKRABber"),
        launch.browser = TRUE
        
    )
}