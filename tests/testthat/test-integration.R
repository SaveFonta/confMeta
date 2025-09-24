# ========================================================================
# Fast check to see if integrate_f works properly
# ========================================================================


test_that("integrate_f on an easy function", {
  f <- function(mu) pmax(0, 1 - abs(mu))
  res <- integrate_f(max_iter = 7, f, lower = -1, upper = 1)
  expect_s3_class(res, "integrate")
  expect_equal(res$value, 1, tolerance = 1e-6)
})
