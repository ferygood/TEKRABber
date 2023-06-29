context("DECorrInputs")

data(fetchDataHmChimp)

inputBundle <- DECorrInputs(fetchDataHmChimp)

test_that("DECorrInputs() returns 2 tables", {
  
  # check there are 6 output tables
  expect_equal(length(inputBundle), 2)
  
})
