context("assay_tekcorrset")

data(speciesCorr)
hmGeneCorrInput <- assay_tekcorrset(speciesCorr, "gene", "human")
hmTECorrInput <- assay_tekcorrset(speciesCorr, "te", "human")
chimpGeneCorrInput <- assay_tekcorrset(speciesCorr, "gene", "chimpanzee")
chimpTECorrInput <- assay_tekcorrset(speciesCorr, "te", "chimpanzee")

data(ctCorr)
geneConCorrInput <- assay_tekcorrset(ctCorr, "gene", "control")
teConCorrInput <- assay_tekcorrset(ctCorr, "te", "control")
geneTreatCorrInput <- assay_tekcorrset(ctCorr, "gene", "treatment")
teTreatCorrInput <- assay_tekcorrset(ctCorr, "te", "treatment")


test_that("assay_tekcorrset() return correct dataframes", {
    
    # check assay_tekcorrset can return the correct dataframes
    expect_vector(hmGeneCorrInput, size = 50)
    expect_vector(hmTECorrInput, size = 50)
    expect_vector(chimpGeneCorrInput, size = 50)
    expect_vector(chimpTECorrInput, size = 50)
    
    expect_vector(geneConCorrInput, size = 10)
    expect_vector(teConCorrInput, size = 10)
    expect_vector(geneTreatCorrInput, size = 10)
    expect_vector(teTreatCorrInput, size = 10)
    
})
