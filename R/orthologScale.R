#' Get orthology information from Ensembl
#' @description Get orthology information from Ensembl using biomaRt and
#' calculate scaling factor between two species using the confidence of
#' orthology score and expression counts.
#' @usage orthologScale(speciesRef, speciesCompare, geneCountRef, 
#' geneCountCompare)
#' @param speciesRef The scientific name for your reference species. 
#' i.e., hsapiens
#' @param speciesCompare The scientific name for your species to compare. 
#' i.e., ptroglodytes
#' @param geneCountRef Gene count from your reference species. First column 
#' should be Ensmebl gene ID
#' @param geneCountCompare Gene count from the species you want to compare. 
#' First column should also be Ensembl gene ID
#' @return There are two outputs:(1) orthologTable: orthology information 
#' from BioMart (2) scale_factor: for normalizing expression counts
#' @export
#' @examples
#' data(speciesCounts)
#' hmGene <- speciesCounts$hmGene
#' chimpGene <- speciesCounts$chimpGene
#' 
#' ## For demonstration, here we only select 1000 rows to save time
#' set.seed(1234)
#' hmGeneSample <- hmGene[sample(nrow(hmGene), 1000), ]
#' chimpGeneSample <- chimpGene[sample(nrow(chimpGene), 1000), ]
#' 
#' fetchData <- orthologScale(
#'     speciesRef = "hsapiens",
#'     speciesCompare = "ptroglodytes",
#'     geneCountRef = hmGeneSample,
#'     geneCountCompare = chimpGeneSample
#' )
orthologScale <- function(
    speciesRef, speciesCompare, 
    geneCountRef, geneCountCompare) {
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
    orthologTable$refLength <- abs(
        orthologTable[["refEnd"]] - orthologTable[["refStart"]])
    
    orthologTable$compareLength <- abs(
        orthologTable[["compareEnd"]] - orthologTable[["compareStart"]])
    
    colnames(geneCountRef)[1] <- "refEnsemblID"
    idx_ref <- seq_len(ncol(geneCountRef))[-1]
    sub_ref <- geneCountRef[, idx_ref]
    geneCountRef$refMean <- rowMeans(sub_ref)
    
    colnames(geneCountCompare)[1] <- "compareEnsemblID"
    idx_compare <- seq_len(ncol(geneCountCompare))[-1]
    sub_compare <- geneCountCompare[, idx_compare]
    geneCountCompare$compareMean <- rowMeans(sub_compare)
    
    ## rearrange data based on confidence of orthology for SCBN estimation
    df <- merge(
        orthologTable,
        geneCountRef[, c(1, ncol(geneCountRef))],
        by = c("refEnsemblID", "refEnsemblID")
    )
    
    df <- merge(
        df,
        geneCountCompare[, c(1, ncol(geneCountCompare))],
        by = c("compareEnsemblID", "compareEnsemblID")
    )
    
    confidence_count <- nrow(df[df[["orthologyConfidence"]]==1, ])
    
    df.scbn <- df[, c("refLength", "refMean", "compareLength", "compareMean")]
    
    ## run scbn to obtain scaling factor
    factor <- SCBN::SCBN(
        orth_gene = df.scbn, 
        hkind = seq_len(confidence_count), 
        a = 0.05)
    
    scaleFactor <- factor$scbn_val
    
    result <- list(
        "orthologTable" = orthologTable,
        "scaleFactor" = scaleFactor
    )
    result
}