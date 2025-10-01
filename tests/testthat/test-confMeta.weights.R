# ========================================================================
# Tests for weights (w) argument with p_edgington_w
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
      fun = p_edgington_w,   # <-- use p_edgington_w directly
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





################################################################################
test_that("Weighted vs unweighted", {
  est <- c(-0.5, -0.1, 0.2)
  se  <- c(0.2, 0.25, 0.3)
  w   <- c(2, 1, 1)
  
  obj_unw <- confMeta(
    estimates = est,
    SEs = se,
    conf_level = 0.95,
    fun = p_edgington_w,
    fun_name = "Edgington unweighted"
  )
  
  obj_w <- confMeta(
    estimates = est,
    SEs = se,
    conf_level = 0.95,
    fun = p_edgington_w,
    fun_name = "Edgington weighted",
    w = w
  )
  
  # what MUST not change
  expect_equal(obj_unw$estimates,       obj_w$estimates)
  expect_equal(obj_unw$SEs,             obj_w$SEs)
  expect_equal(obj_unw$study_names,     obj_w$study_names)
  expect_equal(obj_unw$conf_level,      obj_w$conf_level)
  expect_equal(obj_unw$individual_cis,  obj_w$individual_cis)
  expect_equal(obj_unw$comparison_cis,  obj_w$comparison_cis)
  expect_equal(obj_unw$comparison_p_0,  obj_w$comparison_p_0)
  
  # what MUST change
  expect_false(isTRUE(all.equal(obj_unw$w,          obj_w$w)))
  expect_false(isTRUE(all.equal(obj_unw$joint_cis,  obj_w$joint_cis)))
  expect_false(isTRUE(all.equal(obj_unw$gamma,      obj_w$gamma)))
  expect_false(isTRUE(all.equal(obj_unw$p_max,      obj_w$p_max)))
  expect_false(isTRUE(all.equal(obj_unw$p_0,        obj_w$p_0)))
  expect_false(isTRUE(all.equal(obj_unw$aucc,       obj_w$aucc)))
  expect_false(isTRUE(all.equal(obj_unw$aucc_ratio, obj_w$aucc_ratio)))
})


##################################################################à
test_that("Weighted with weights=1 vs unweighted", {
  est <- c(-0.5, -0.1, 0.2)
  se  <- c(0.2, 0.25, 0.3)
  w1  <- rep(1, length(est))
  
  obj_unw <- confMeta(
    estimates = est,
    SEs = se,
    conf_level = 0.95,
    fun = p_edgington_w,
    fun_name = "Edgington unweighted"
  )
  
  obj_w1 <- confMeta(
    estimates = est,
    SEs = se,
    conf_level = 0.95,
    fun = p_edgington_w,
    fun_name = "Edgington weighted (all ones)",
    w = w1
  )
  
  # What changes
  ignore <- c("w", "fun_name")
  
  #The rest MUST be the same
  for (nm in setdiff(names(obj_unw), ignore)) {
    expect_equal(obj_unw[[nm]], obj_w1[[nm]])
  }
})







#######################


test_that("The order doesnt change the output", {
  est <- c(-0.5, -0.1, 0.2)
  se  <- c(0.2, 0.25, 0.3)
  w   <- c(2, 1, 3)
  
  obj1 <- confMeta(
    estimates = est,
    SEs = se,
    w = w,
    conf_level = 0.95,
    fun = p_edgington_w
  )
  
  o <- c(3,1,2)  # reorder in a different way
  obj2 <- confMeta(
    estimates = est[o],
    SEs = se[o],
    w = w[o],
    conf_level = 0.95,
    fun = p_edgington_w
  )
  
  #expect equal
  expect_equal(obj1$joint_cis,     obj2$joint_cis)
  expect_equal(obj1$gamma,         obj2$gamma)
  expect_equal(obj1$p_max,         obj2$p_max)
  expect_equal(obj1$p_0,           obj2$p_0)
  expect_equal(obj1$aucc,          obj2$aucc)
  expect_equal(obj1$aucc_ratio,    obj2$aucc_ratio)
  expect_equal(obj1$comparison_cis,obj2$comparison_cis)
  expect_equal(obj1$comparison_p_0,obj2$comparison_p_0)

})
