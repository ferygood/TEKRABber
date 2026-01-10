#' Prepare a table from two species RepeatMakser track from UCSC genome Table
#' @description create a table to the rmsk argument in orthologScale(). 
#' Before version 1.8, TEKRABber requires user to prepare this table by themselves and 
#' this function can help user automatically get the RepeatMasker table from 
#' UCSC. The arguments required are the abbreviation of the version of 
#' reference (case-sensitive). For example, "hg38" for human. 
#' Note: currently only 91 genomes provided. Check if the reference exists with 
#' GenomeInfoDb::registered_UCSC_genomes().
#' 
#' @param refSpecies the version of reference species, i.e. hg38
#' @param compareSpecies the version of compared species, i.e. panTro6
#'
#' @export
#' @return Dataframe with four columns: repName, repClass, rLen and cLen
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom dplyr mutate select group_by summarise
#' @importFrom magrittr %>%
#' @examples
#' df_rmsk <- prepareRMSK(refSpecies = "hg38", compareSpecies = "panTro6")
#' 
prepareRMSK <- function(refSpecies, compareSpecies) {

    fetch_rmsk <- function(species) {
        ah <- AnnotationHub::AnnotationHub()
        q <- AnnotationHub::query(
            ah,
            c("RepeatMasker", "UCSC", species)
        )

        if (length(q) == 0) {
            stop("No RepeatMasker data found for species: ", species)
        }

        rmsk <- q[[1]]
        as.data.frame(rmsk)
    }

    summarise_rmsk <- function(df, len_name) {
        df %>%
            mutate(.len = abs(repEnd - repStart)) %>%
            select(repName, repClass, .len) %>%
            group_by(repName, repClass) %>%
            summarise(
                !!len_name := mean(.len),
                .groups = "drop"
            )
    }

    ref_df <- fetch_rmsk(refSpecies)
    cmp_df <- fetch_rmsk(compareSpecies)

    ref_tbl <- summarise_rmsk(ref_df, "rLen")
    cmp_tbl <- summarise_rmsk(cmp_df, "cLen")

    merge(
        ref_tbl,
        cmp_tbl,
        by = c("repName", "repClass")
    )
}
