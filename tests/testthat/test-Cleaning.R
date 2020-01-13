test_that("Beta-matrix has been input in the correct format in clean_beta", {
  # Check that the row names are valid
  data("GSE121961", envir = environment())
  rownames(GSE121961) <- NULL
  expect_error(
    clean_beta(beta = GSE121961),
    "Please make sure your beta-matrix has row names corresponding to CpGs"
  )
})
