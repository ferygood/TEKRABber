context("DECorrInput")

data(fetchDataHmChimp, speciesCounts)

inputBundle <- DECorrInputs(
  orthologTable = fetchDataHmChimp$orthologTable,
  scaleFactor = fetchDataHmChimp$scaleFactor,
  geneCountRef = speciesCounts$hmGene,
  geneCountCompare = speciesCounts$chimpGene,
  teCountRef = speciesCounts$hmTE,
  teCountCompare = speciesCounts$chimpTE
)

test_that("DECorrInputs() returns 6 tables", {
  
  # check there are 6 output tables
  expect_equal(length(inputBundle), 6)
  
})