# ========================================================================
# Tests for weights (w) argument with p_edgington_w2
# ========================================================================

test_that("weights argument is validated correctly", {
  
  n <- 5
  estimates <- rnorm(n)
  SEs <- rgamma(n, 5, 5)
  conf_level <- 0.95
  
  cl1 <- quote({
    confMeta(
      estimates = estimates,
      SEs = SEs,
      conf_level = conf_level,
      fun = p_edgington_w2,   # <-- use p_edgington_w2 directly
      w = w
    )
  })
  
  # wrong length
  w <- as.numeric(1:4)
  expect_error(eval(cl1))
  
  # contains NA
  w <- c(1, NA, 1, 1, 1)
  expect_error(eval(cl1))
  
  # contains NaN
  w <- c(1, 1, as.numeric(NaN), 1, 1)
  expect_error(eval(cl1))
  
  # contains negative value
  w <- c(-1, 1, 1, 1, 1)
  expect_error(eval(cl1))
  
  # all zero weights (not allowed)
  w <- rep(0, n)
  expect_error(eval(cl1))
  
  # valid integer weights
  w <- as.numeric(1:5)
  expect_s3_class(eval(cl1), "confMeta")
  
  # valid numeric weights (non-integers)
  w <- runif(n, min = 0.1, max = 2)
  expect_s3_class(eval(cl1), "confMeta")
  
  # very large weights
  w <- rep(1e6, n)
  expect_s3_class(eval(cl1), "confMeta")
  
  # very small positive weights
  w <- rep(1e-6, n)
  expect_s3_class(eval(cl1), "confMeta")
})

