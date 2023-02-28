orthologScale_beta <- function(
        speciesRef, speciesCompare, 
        geneCountRef, geneCountCompare,
        teCountRef, teCountCompare, 
        rmsk) {
    
    # Part1: normalize orthologous genes between species
    ## 1-1: Get ortholog table using biomaRt
    geneRef <- paste0(speciesRef, "_gene_ensembl")
    geneCompare <- paste0(speciesCompare, "_gene_ensembl")
    orthologyRef <- paste0(speciesRef, "_homolog_orthology_confidence")
    
    ensemblRef <- biomaRt::useEnsembl("ensembl", dataset = geneRef, version = "105")
    ensemblCompare <- biomaRt::useEnsembl("ensembl", dataset = geneCompare, version = "105")
    
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
    
    ## 1-2: Estimate scaling factor based on mean expression of genes
    ## and orthologTable from 1-1.
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
    
    ## 1-3 rearrange data based on confidence of orthology for SCBN estimation
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
    
    # arrange based on orthologyConfidence
    df <- arrange(df, desc(orthologyConfidence))
    
    confidence_count <- nrow(df[df[["orthologyConfidence"]]==1, ])
    
    df.scbn <- df[, c("refLength", "refMean", "compareLength", "compareMean")]
    
    ## run scbn to obtain scaling factor
    factor <- SCBN::SCBN(
        orth_gene = df.scbn, 
        hkind = seq_len(confidence_count), 
        a = 0.05)
    
    c_ortholog <- factor$scbn_val  #scaling factor
    
    ## 1-4 normalize data
    N1 <- sum(df.scbn$refMean)
    N2 <- sum(df.scbn$compareMean)
    
    # remove mean from previous step
    geneCountCompare <- geneCountCompare[, -ncol(geneCountCompare)]
    
    geneCountCompare <- merge(
        geneCountCompare,
        df[,c("compareEnsemblID", "refLength", "compareLength")],
        by = c("compareEnsemblID", "compareEnsemblID")
    )
    
    normalized_geneCountCompare <- geneCountCompare %>%
        mutate(scale = (refLength * N1 * c_ortholog)/(compareLength * N2)) %>%
        mutate(across(2:(ncol(geneCountCompare)-2), ~ .x * scale))
    
    normalized_geneCountCompare <- 
        normalized_geneCountCompare[,c(1:(ncol(normalized_geneCountCompare)-3))]
    
    geneCountRef <- geneCountRef[geneCountRef$refEnsemblID %in% df$refEnsemblID, 
                                 -ncol(geneCountRef)]
    
    
    ## Part2: normalize TE between species
    # 2-1 Calculate mean expression and mean gene length into a dataframe
    teCountRef$refMean <- rowMeans(teCountRef[, -1])
    teCountCompare$compareMean <- rowMeans(teCountCompare[, -1])
    
    # rename column name for merging dataframe
    colnames(teCountRef)[1] <- "teName"
    colnames(teCountCompare)[1] <- "teName"
    colnames(rmsk) <- c("teName", "teClass", "refLen", "compareLen")
    
    df.TE <- merge(rmsk, teCountRef[,c("teName", "refMean")],
                   by = c("teName", "teName"))
    df.TE <- merge(df.TE, teCountCompare[,c("teName", "compareMean")],
                   by = c("teName", "teName"))
    
    # rearrange order by select_class
    select_class <- c("LTR", "LINE", "SINE", "DNA")
    df.TE.conf <- df.TE[df.TE$teClass %in% select_class, ]
    df.TE.notconf <- df.TE[!(df.TE$teClass %in% select_class), ]
    
    df.TE.combine <- rbind(df.TE.conf, df.TE.notconf)
    
    # 2-2 Calcualte scaling factor 
    te.factor <- SCBN::SCBN(
        orth_gene = 
            df.TE.combine[,c("refLen", "refMean", "compareLen", "compareMean")], 
        hkind = seq_len(nrow(df.TE.conf)), 
        a = 0.05)
    
    c_te <- te.factor$scbn_val  #scaling factor
    
    print(c_te)
    # 2-3 Normalize data
    T1 <- sum(df.TE.combine$refMean)
    T2 <- sum(df.TE.combine$refCompare)
    
    teCountCompare <- teCountCompare[, -ncol(teCountCompare)]
    
    teCountCompare <- merge(
        teCountCompare,
        df.TE.combine[,c("teName", "refLen", "compareLen")],
        by = c("teName", "teName")
    )
    
    normalized_teCountCompare <- teCountCompare %>%
        mutate(scale = (refLen * T1 * c_te)/(compareLen * T2)) %>%
        mutate(across(2:(ncol(teCountCompare)-2), ~ .x * scale))
    
    normalized_teCountCompare <- 
        normalized_teCountCompare[,c(1:(ncol(normalized_teCountCompare)-3))]
    
    teCountRef <- teCountRef[teCountRef$teName %in% df.TE.combine$teName, 
                                 -ncol(teCountRef)]
    
    result <- list(
        "orthologTable" = orthologTable,
        "c_ortholog" = c_ortholog,
        "geneRef" = geneCountRef,
        "geneCompare" = normalized_geneCountCompare,
        "c_te" = c_te,
        "teRef" = teCountRef,
        "teCompare" = normalized_teCountCompare
    )

    result
    
    
}
