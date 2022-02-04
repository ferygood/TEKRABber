context("orthologScale")

data(speciesCounts)
hmGene <- speciesCounts$hmGene
hmTE <- speciesCounts$hmTE
chimpGene <- speciesCounts$chimpGene
chimpTE <- speciesCounts$chimpTE

fetchData <- orthologScale(
  geneCountRef = hmGene,
  geneCountCompare = chimpGene,
  speciesRef = "hsapiens",
  speciesCompare = "ptroglodytes"
)

col_names <- c("refGene", "chromosome", "refEnsemblID", "refStart", "refEnd",
              "compareGene", "compareEnsemblID", "compareStart", 
              "compareEnd", "orthologyConfidence", "refLength", "compareLength")


test_that("orthologScale() gives correct table and one scaling factor", {
  
  # check column names
  expect_named(fetchData$orthologTable, col_names)
  
  # check there are two outputs
  expect_equal(length(fetchData), 2)

})