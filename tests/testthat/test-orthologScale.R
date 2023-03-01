context("orthologScale")

data(speciesCounts)
data(hg38_panTro6_rmsk)
hmGene <- speciesCounts$hmGene
chimpGene <- speciesCounts$chimpGene
hmTE <- speciesCounts$hmTE
chimpTE <- speciesCounts$chimpTE

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
        geneCountCompare = chimpGeneSample,
        teCountRef = hmTE,
        teCountCompare = chimpTE,
        rmsk = hg38_panTro6_rmsk
    )
    
    expect_equal(length(fetchData), 7)
    expect_s3_class(fetchData$orthologTable, "data.frame")
    expect_type(fetchData$c_ortholog, "double")
    expect_type(fetchData$c_te, "double")
    expect_named(fetchData$orthologTable, col)
})
