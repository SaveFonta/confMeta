# ========================================================================
# Regression test for p_edgington_w with known values
# ========================================================================

test_that("p_edgington_w reproduces known CI for fixed inputs", {
  
  estimates <- c( 0.69314718, -0.22292056, -0.53382746,
                  -0.78642817,  1.38629436, -0.33551242, -0.09531018)
  
  SEs <- c(1.1401754, 0.2526471, 0.1446606,
           0.4188102, 0.9279607, 0.3185523, 0.5871586)
  
  w <- c(0.04373144, 0.19735637, 0.34467926,
          0.11905516, 0.05373236, 0.15652540, 0.08492001)
  
  # compute CI using your Edgington weighted function
  cm <- confMeta(
    estimates = estimates,
    SEs = SEs,
    conf_level = 0.95,
    fun = p_edgington_w,
    w = w,
    input_p = "greater",
    alternative = "one.sided",
    heterogeneity = "none"
  )
  
  ci <- cm$joint_cis
  
  # expected CI from reference run
  expected <- matrix(c(-0.6085842, -0.1325605), nrow = 1,
                     dimnames = list(NULL, c("lower", "upper")))
  
  # check that CI matches within numerical tolerance
  expect_equal(ci, expected, tolerance = 1e-6)
})



# ========================================================================
# Other tests more standard just to check if the w is valid as an argument
# ========================================================================

test_that("p_edgington_w returns valid p-values", {
  
  n <- 5
  estimates <- rnorm(n)
  SEs <- rgamma(n, 5, 5)
  w <- rep(1, n)
  
  # with default weights (all 1)
  pval <- p_edgington_w(estimates, SEs, mu = 0, w = w)
  expect_true(is.numeric(pval))
  expect_true(all(pval >= 0 & pval <= 1))
  
  # with custom weights
  w <- as.numeric(1:n)
  pval_w <- p_edgington_w(estimates, SEs, mu = 0, w = w)
  expect_true(is.numeric(pval_w))
  expect_true(all(pval_w >= 0 & pval_w <= 1))
  
  # p-value should be finite
  expect_true(all(is.finite(pval_w)))
})

# ------------------------------------------------------------------------
# Comparison to unweighted Edgington (p_edgington)
# ------------------------------------------------------------------------

test_that("p_edgington_w equals p_edgington when all weights = 1", {
  
  n <- 5
  estimates <- rnorm(n)
  SEs <- rgamma(n, 5, 5)
  
  p1 <- p_edgington(estimates, SEs, mu = 0)
  p2 <- p_edgington_w(estimates, SEs, mu = 0, w = rep(1, n))
  
  # they should be equal (or numerically very close)
  expect_equal(p1, p2, tolerance = 1e-12)
})

# ------------------------------------------------------------------------
# Edge cases for p_edgington_w
# ------------------------------------------------------------------------

test_that("edge cases for weights in p_edgington_w", {
  
  estimates <- rnorm(3)
  SEs <- runif(3, min = 0.5, max = 1.5)
  
  
  
  # mismatched length between estimates and weights
  expect_error(p_edgington_w(estimates, SEs, mu = 0, w = c(1, 2)))
  
  # NA in weights
  expect_error(p_edgington_w(estimates, SEs, mu = 0, w = c(1, NA, 1)))
  
  # negative weight
  expect_error(p_edgington_w(estimates, SEs, mu = 0, w = c(-1, 2, 3)))
  
  # zero weights (not allowed)
  expect_error(p_edgington_w(estimates, SEs, mu = 0, w = c(0, 1, 2)))
  
  # very large weights
  p_large <- p_edgington_w(estimates, SEs, mu = 0, w = c(1e6, 1e6, 1e6))
  expect_true(is.numeric(p_large))
  expect_true(all(p_large >= 0 & p_large <= 1))
  
  # very small positive weights
  p_small <- p_edgington_w(estimates, SEs, mu = 0, w = c(1e-6, 1e-6, 1e-6))
  expect_true(is.numeric(p_small))
  expect_true(all(p_small >= 0 & p_small <= 1))
})
