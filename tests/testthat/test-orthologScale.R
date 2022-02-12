context("orthologScale")

data(speciesCounts)
hmGene <- speciesCounts$hmGene
chimpGene <- speciesCounts$chimpGene

set.seed(1234)
hmGeneSample <- hmGene[sample(nrow(hmGene), 1000), ]
chimpGeneSample <- chimpGene[sample(nrow(chimpGene), 1000), ]

# expect column names
col <- c("refGene", "chromosome", "refEnsemblID", "refStart", "refEnd", 
         "compareGene", "compareEnsemblID", "compareStart", "compareEnd",
         "orthologyConfidence", "refLength", "compareLength")

test_that("orthologScale returns orthology information and scaling factor", {
    
    fetchData <- orthologScale(
        speciesRef = "hsapiens",
        speciesCompare = "ptroglodytes",
        geneCountRef = hmGeneSample,
        geneCountCompare = chimpGeneSample
    )
    
    expect_equal(length(fetchData), 2)
    expect_s3_class(fetchData$orthologTable, "data.frame")
    expect_type(fetchData$scaleFactor, "double")
    expect_named(fetchData$orthologTable, col)
})
