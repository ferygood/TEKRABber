#' Get orthology information from Ensembl
#' @description Get orthology information from Ensembl using biomaRt and
#' calculate scaling factor between two species using the confidence of
#' orthology score and expression counts.
#' @usage orthologScale(speciesRef, speciesCompare, geneCountRef, geneCountCompare)
#' @param speciesRef The scientific name for your reference species. i.e. hsapiens
#' @param speciesCompare The scientific name for your species to compare. i.e. ptroglodytes
#' @param geneCountRef Gene count from your reference species. First column should be Ensmebl gene ID
#' @param geneCountCompare Gene count from the species you want to compare. First column should also be Ensembl gene ID
#' @return There are two outputs:(1) orthologTable: orthology information from BioMart (2) scale_factor: for normalizing expression counts
#' @importFrom dplyr filter
#' @export
#' @examples
#' data(speciesCounts)
#' hmGene <- speciesCounts$hmGene
#' chimpGene <- speciesCounts$chimpGene
#' #fetchData <- orthologScale(
#' #  speciesRef = "hsapiens",
#' #  speciesCompare = "ptroglodytes",
#' #  geneCountRef = hmGene,
#' #  geneCountCompare = chimpGene
#' #)
orthologScale <- function(speciesRef, speciesCompare, geneCountRef, geneCountCompare) {
    ## Part1: Get ortholog table using biomaRt
    geneRef <- paste0(speciesRef, "_gene_ensembl")
    geneCompare <- paste0(speciesCompare, "_gene_ensembl")
    orthologyRef <- paste0(speciesRef, "_homolog_orthology_confidence")
    
    ensemblRef <- biomaRt::useEnsembl("ensembl", dataset = geneRef)
    ensemblCompare <- biomaRt::useEnsembl("ensembl", dataset = geneCompare)
    
    orthologTable <- biomaRt::getLDS(
        attributes = c(
            "external_gene_name",
            "chromosome_name",
            "ensembl_gene_id",
            "start_position",
            "end_position"
        ),
        mart = ensemblRef,
        attributesL = c(
            "external_gene_name",
            "ensembl_gene_id",
            "start_position",
            "end_position",
            orthologyRef
        ),
        martL = ensemblCompare
    )
    
    colnames(orthologTable) <- c(
        "refGene",
        "chromosome",
        "refEnsemblID",
        "refStart",
        "refEnd",
        "compareGene",
        "compareEnsemblID",
        "compareStart",
        "compareEnd",
        "orthologyConfidence"
    )
    
    ## Part2: Estimate scaling factor based on mean expression of genes
    ## and orthologTable from Part1.
    orthologTable <- orthologTable %>%
        mutate(refLength = refEnd - refStart) %>%
        mutate(compareLength = compareEnd - compareStart)
    
    colnames(geneCountRef)[1] <- "refEnsemblID"
    geneCountRef <- geneCountRef %>%
        mutate(refMean = rowMeans(select(., 2:ncol(geneCountRef))))
    
    colnames(geneCountCompare)[1] <- "compareEnsemblID"
    geneCountCompare <- geneCountCompare %>%
        mutate(compareMean = rowMeans(select(., 2:ncol(geneCountCompare))))
    
    ## rearrange data based on confidence of orthology for SCBN estimation
    df <- inner_join(
        orthologTable,
        geneCountRef[, c(1, ncol(geneCountRef))],
        by = "refEnsemblID"
    )
    
    df <- inner_join(
        df,
        geneCountCompare[, c(1, ncol(geneCountCompare))],
        by = "compareEnsemblID"
    )
    
    confidence_count <- df %>%
        filter(orthologyConfidence == 1) %>%
        nrow()
    
    df.scbn <- df %>% select(c(11, 13, 12, 14))
    
    ## run scbn to obtain scaling factor
    factor <- SCBN::SCBN(orth_gene = df.scbn, hkind = 1:confidence_count, a = 0.05)
    scaleFactor <- factor$scbn_val
    
    result <- list(
        "orthologTable" = orthologTable,
        "scaleFactor" = scaleFactor
    )
    result
}