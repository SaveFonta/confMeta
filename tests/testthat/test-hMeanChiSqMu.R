test_that("Results of hMeanChiSqMu are consistent.", {
  # Set inputs
  set.seed(42)
  n <- 15
  thetahat <- rnorm(n)
  se <- rgamma(n, 5, 5)
  mu <- seq(
    min(thetahat) - 0.5 * max(se),
    max(thetahat) + 0.5 * max(se),
    length.out = 100
  )
  phi <- estimatePhi(thetahat = thetahat, se = se)
  tau2 <- estimateTau2(thetahat = thetahat, se = se)
  grid <- expand.grid(
    heterogeneity = c("none", "additive", "multiplicative"),
    distr = c("chisq", "f"),
    alternative = c("none", "less", "greater", "two.sided"),
    stringsAsFactors = FALSE
  )
  # Run new function
  res <- lapply(
    seq_len(nrow(grid)),
    function(x) {
      het <- grid$heterogeneity[x]
      phi <- if (het == "multiplicative") phi else NULL
      tau2 <- if (het == "additive") tau2 else NULL
      distr <- grid$distr[x]
      alternative <- grid$alternative[x]
      hMeanChiSqMu(
        thetahat = thetahat,
        se = se,
        mu = mu,
        phi = phi,
        tau2 = tau2,
        heterogeneity = het,
        alternative = alternative,
        check_inputs = TRUE,
        distr = distr
      )
    }
  )
  
  # Get the old function, vectorise it, and run the same inputs
  old_fun <- get_old_FUN(
    path = "https://raw.githubusercontent.com/felix-hof/hMean/main/R/hMeanChiSqMu.R",
    fun_name = "hMeanChiSqMu"
  )
  
  old_res <- lapply(
    seq_len(nrow(grid)),
    function(x) {
      het <- grid$heterogeneity[x]
      phi <- if (het == "multiplicative") phi else NULL
      tau2 <- if (het == "additive") tau2 else NULL
      distr <- grid$distr[x]
      alternative <- grid$alternative[x]
      old_fun(
        thetahat = thetahat,
        se = se,
        mu = mu,
        phi = phi,
        tau2 = tau2,
        heterogeneity = het,
        alternative = alternative,
        check_inputs = TRUE,
        distr = distr
      )
    }
  )
  
  # compare results
  expect_equal(res, old_res)
})
