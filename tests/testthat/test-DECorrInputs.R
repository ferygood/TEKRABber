context("DECorrInput")

data(fetchDataHmChimp)

inputBundle <- DECorrInputs(fetchDataHmChimp)

test_that("DECorrInputs() returns 6 tables", {
  
  # check there are 6 output tables
  expect_equal(length(inputBundle), 6)
  
})
