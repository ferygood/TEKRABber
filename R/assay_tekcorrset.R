#' Access gene and transposable elements expression data
#' @description Access expression data from the TekCorrSet class for running vignette
#' demonstration.
#' @param tecorrset TekCorrSet object 
#' @param expType Indicate which data you want to access. It should be "gene" or "te".
#' @param sample The species name or experimental design. It should be "human" and 
#' "chimpanzee" when you are running comparing species design. "control" and "treatment"
#'  are for running same species design.
#' @return a dataframe contains expression genes or transposable element.
#' @export
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
            sprintf("Wrong argument. It should be %s or %s.", sepc[1,1], spec[2,1])
        )
    } else if (expType == "te") {
        if (sample == spec[1,1]){
            assay(tecorrset@teSE[, tecorrset@teSE$sample==spec[1,1]])
        } else if (sample == spec[2,1]) {
            assay(tecorrset@teSE[, tecorrset@teSE$sample==spec[2,1]])
        } else (
            sprintf("Wrong argument. It should be %s or %s.", sepc[1,1], spec[2,1])
        )
    } else {
        "Wrong argument. It should be gene or te."
    }
}