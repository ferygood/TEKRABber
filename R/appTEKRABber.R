#' Visualize TEKRABber results with shiny app
#' @description To help user explore their results using TEKRABber, this 
#' function visualizes the results using a self-written shiny app with two 
#' tabs, including the expression and correlation of genes and TEs. To run 
#' it, you need to create four variables and assign them with your DE result, 
#' correlation results and metadata to appDE, appRef, appCompare and 
#' appMeta. Please see the example below for more details.
#' @usage appTEKRABber()
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
#' # create metadata for DE analysis
#' meta <- data.frame(species=c(rep("human", ncol(hmGene) - 1), 
#'     rep("chimpanzee", ncol(chimpGene) - 1))
#' )
#' rownames(meta) <- colnames(inputBundle$geneInputDESeq2)
#' meta$species <- factor(meta$species, levels = c("human", "chimpanzee"))
#' 
#' # DE analysis
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
#' 
#' # Correlation analysis
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
#' # assign results and metadata to appDE, appRef, appCompare, and appMeta
#' appDE <- hmchimpDE
#' appRef <- hmCorrResult
#' appCompare <- chimpCorrResult
#' appMeta <- meta
#' 
#' if (interactive()){
#'     
#'     appTEKRABber()
#' 
#' }
#' 
appTEKRABber <- function() {
    
    shiny::runApp(
        appDir = system.file("shinyGUI", package="TEKRABber"),
        launch.browser = TRUE
        
    )
}