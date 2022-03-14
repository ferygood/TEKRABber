context("DEgeneTE")

data(ctInputDE)
geneInputDE <- ctInputDE$gene
teInputDE <- ctInputDE$te

metaExp <- data.frame(experiment = c(rep("control", 5), rep("treatment", 5)))
rownames(metaExp) <- colnames(geneInputDE)
metaExp$experiment <- factor(metaExp$experiment, 
                             levels = c("control", "treatment"))

resultDE <- DEgeneTE(
    geneTable = geneInputDE,
    teTable = teInputDE,
    metadata = metaExp,
    expDesign = FALSE
)


test_that("DEgeneTE() returns the correct dds and res tables", {
    
    expect_s4_class(resultDE$gene_dds, "DESeqDataSet")
    expect_s4_class(resultDE$gene_res, "DESeqResults")
    expect_s4_class(resultDE$te_dds, "DESeqDataSet")
    expect_s4_class(resultDE$te_res, "DESeqResults")
    
})
