#' Access genes and transposable elements expression data
#' @description a function only used for accessing the expression data from a 
#' TekCorrSet class object to demonstrate examples in vignettes.
#' demonstration.
#' @param tecorrset TekCorrSet object 
#' @param expType Indicate which data you want to access. It should be "gene" 
#' or "te".
#' @param sample The species name or experimental design. It should be "human" 
#' and "chimpanzee" when you are running comparing species design. "control" 
#' and "treatment" are for running same species design.
#' @return a dataframe contains expression genes or transposable elements.
#' @export
#' @importFrom SummarizedExperiment colData assay
#' @examples
#' data(speciesCorr)
#' hmGeneCorrInput <- assay_tekcorrset(speciesCorr, "gene", "human")
#' hmTECorrInput <- assay_tekcorrset(speciesCorr, "te", "human")
#' chimpGeneCorrInput <- assay_tekcorrset(speciesCorr, "gene", "chimpanzee")
#' chimpTECorrInput <- assay_tekcorrset(speciesCorr, "te", "chimpanzee")
assay_tekcorrset <- function(tecorrset, expType, sample) {
    spec <- unique(colData(tecorrset@geneSE))
    if (expType == "gene") {
        if (sample == spec[1,1]){
            assay(tecorrset@geneSE[, tecorrset@geneSE$sample==spec[1,1]])
        } else if (sample == spec[2,1]) {
            assay(tecorrset@geneSE[, tecorrset@geneSE$sample==spec[2,1]])
        } else (
            sprintf(
                "Wrong parameter. It should be %s or %s.", 
                spec[1,1], spec[2,1]
            )
        )
    } else if (expType == "te") {
        if (sample == spec[1,1]){
            assay(tecorrset@teSE[, tecorrset@teSE$sample==spec[1,1]])
        } else if (sample == spec[2,1]) {
            assay(tecorrset@teSE[, tecorrset@teSE$sample==spec[2,1]])
        } else (
            sprintf(
                "Wrong parameter. It should be %s or %s.", 
                spec[1,1], spec[2,1]
            )
        )
    } else {
        "Wrong parameter. It should be gene or te."
    }
}