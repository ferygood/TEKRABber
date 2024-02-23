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
#' @return Dataframe with four columns: repName, repClass, rLen and cLen
#' @export
#' @importFrom rtracklayer browserSession getTable ucscTableQuery
#' @importFrom GenomeInfoDb genome
#' @importFrom dplyr mutate select group_by summarise
#' @importFrom magrittr %>%
#' @examples 
#' df_rmsk <- prepareRMSK(refSpecies = "hg38", compareSpecies = "panTro6") 
#' 
prepareRMSK <- function(refSpecies, compareSpecies){
    # create a session and query repeatmakser track
    # reference species
    refSession <- browserSession("UCSC")
    GenomeInfoDb::genome(refSession) <- refSpecies
    ref.rmsk <- getTable(
        ucscTableQuery(
            refSession, 
            track="RepeatMasker", 
            table="rmsk")
    )
    
    ref.rmsk.tbl <- ref.rmsk %>%
        mutate(rLen = abs(repEnd - repStart)) %>%
        select(c(repName, repClass, rLen)) %>%
        group_by(repName, repClass) %>%
        summarise(rLen = abs(mean(rLen)))
    
    # compare species
    compareSession <- browserSession("UCSC")
    GenomeInfoDb::genome(compareSession) <- compareSpecies
    compare.rmsk <- getTable(
        ucscTableQuery(
            compareSession, 
            track="RepeatMasker", 
            table="rmsk")
    )
    
    compare.rmsk.tbl <- compare.rmsk %>%
        mutate(cLen = abs(repEnd - repStart)) %>%
        select(c(repName, repClass, cLen)) %>%
        group_by(repName, repClass) %>%
        summarise(rLen = abs(mean(cLen)))
    
    # merge table
    merged_df <- merge(ref.rmsk.tbl, compare.rmsk.tbl, by=c("repName", "repClass"))
    colnames(merged_df)[c(3,4)] <- c("rLen", "cLen")
    
    merged_df
}
