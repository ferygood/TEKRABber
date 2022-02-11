context("orthologScale")

data(fetchDataHmChimp)

# expect column names
col <- c("refGene", "chromosome", "refEnsemblID", "refStart", "refEnd", 
         "compareGene", "compareEnsemblID", "compareStart", "compareEnd",
         "orthologyConfidence", "refLength", "compareLength")


test_that("orthologScale returns orthology information and scaling factor", {
    
    fetchData <- fetchDataHmChimp
    
    expect_equal(length(fetchData), 2)
    expect_s3_class(fetchData$orthologTable, "data.frame")
    expect_type(fetchData$scaleFactor, "double")
    expect_named(fetchData$orthologTable, col)
})
