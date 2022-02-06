context("corrOrthologTE")

data(ctCorr)
geneConCorrInput <- assay_tekcorrset(ctCorr, "gene", "control")
teConCorrInput <- assay_tekcorrset(ctCorr, "te", "control")

controlCorr <- corrOrthologTE(
    geneInput = geneConCorrInput,
    teInput = teConCorrInput,
    corrMethod = "pearson",
    padjMethod = "fdr"
)

test_that("corrOrthologTE return the correct correlation result", {
    local_edition(3)
    expect_snapshot(controlCorr)
    
})
