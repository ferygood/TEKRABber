context("corrOrthologTE")

# load built-in data
data(ctInputDE)
geneInputDE <- ctInputDE$gene
teInputDE <- ctInputDE$te

metaExp <- data.frame(experiment = c(rep("control", 5), rep("treatment", 5)))
rownames(metaExp) <- colnames(geneInputDE)
metaExp$experiment <- factor(
    metaExp$experiment, 
    levels = c("control", "treatment")
)

resultDE <- DEgeneTE(
    geneTable = geneInputDE,
    teTable = teInputDE,
    metadata = metaExp,
    expDesign = FALSE
)

controlCorr <- corrOrthologTE(
    geneInput = resultDE$geneCorrInputRef[c(1:10),],
    teInput = resultDE$teCorrInputRef[c(1:10),],
    corrMethod = "pearson",
    padjMethod = "fdr"
)

test_that("corrOrthologTE return the correct correlation result", {
    local_edition(3)
    expect_snapshot(controlCorr)
    
})
