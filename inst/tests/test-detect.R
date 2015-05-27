context("Anomalous time series tests")

test_that("Time series is constructed correctly, features are extracted properly and that anomalous method returns expected number of outliers.", {
  z <- ts(matrix(rnorm(3000),ncol=100),freq=4)
  expect_that(ncol(z) == 100, is_true())
  y <- tsmeasures(z)
  expect_that(ncol(y) == 16, is_true())
  expect_that(length(anomaly(y, n=2, method="ahull", plot=FALSE)) == 2, is_true()) 
  expect_that(length(anomaly(y, n=2, method="hdr", plot=FALSE)) == 2, is_true()) 
})
