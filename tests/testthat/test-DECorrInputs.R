context("DECorrInput")

data(fetchDataHmChimp, speciesCounts)
fetchData <- fetchDataHmChimp

inputBundle <- DECorrInputs(
  orthologTable = fetchData$orthologTable,
  scaleFactor = fetchData$scaleFactor,
  geneCountRef = hmGene,
  geneCountCompare = chimpGene,
  teCountRef = hmTE,
  teCountCompare = chimpTE
)

test_that("DECorrInputs() returns 6 tables", {
  
  # check there are 6 output tables
  expect_equal(length(inputBundle), 6)
  
})