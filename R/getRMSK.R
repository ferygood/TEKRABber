#' Get reference and compare species repeatmakser from UCSC genome Table
#'
#' @param refSpecies the version of reference, i.e. hg38
#' @param compareSpecies the version of comparison, i.e. panTro6
#'
#' @return Dataframe with four columns, repName, repClass, rLen and cLen
#' @export
#'
#' @examples 
#' 
#' df_rmsk <- getRMSK(refSpecies = "hg38", compareSpecies = "panTro6") 
#' 
getRMSK <- function(refSpecies, compareSpecies){
    # create a session and query repeatmakser track
    # reference species
    refSession <- rtracklayer::browserSession("UCSC")
    rtracklayer::genome(refSession) <- refSpecies
    ref.rmsk <- rtracklayer::getTable(
        rtracklayer::ucscTableQuery(
            refSession, 
            track="RepeatMasker", 
            table="rmsk")
    )
    
    ref.rmsk.tbl <- ref.rmsk %>%
        dplyr::mutate(rLen = abs(repEnd - repStart)) %>%
        dplyr::select(c(repName, repClass, rLen)) %>%
        dplyr::group_by(repName, repClass) %>%
        dplyr::summarise(rLen = abs(mean(rLen)))
    
    # compare species
    compareSession <- rtracklayer::browserSession("UCSC")
    rtracklayer::genome(compareSession) <- comapreSpecies
    compare.rmsk <- rtracklayer::getTable(
        rtracklayer::ucscTableQuery(
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
    merge_df <- merge(ref.rmsk.tbl, compare.rmsk.tbl, by=c("repName", "repClass"))
    colnames(merge_df)[3,4] <- c("rLen", "cLen")
    
    merged_df
}
