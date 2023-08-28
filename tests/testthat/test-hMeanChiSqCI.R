# source globals
source("globals.R")

test_that("Results of hMeanChiSqCI are consistent.", {
  # Set inputs
  set.seed(42)
  n <- 15
  thetahat <- rnorm(n)
  se <- rgamma(n, 5, 5)
  phi <- estimatePhi(thetahat = thetahat, se = se)
  tau2 <- estimateTau2(thetahat = thetahat, se = se)
  wGamma <- rep(1, length(thetahat) - 1)
  level <- 0.99
  check_inputs <- TRUE

  grid <- expand.grid(
    alternative = c("none", "less", "greater", "two.sided"),
    pValueFUN = c("hMeanChiSqMu", "kTRMu", "pPearsonMu"),
    heterogeneity = c("none", "additive", "multiplicative"),
    stringsAsFactors = FALSE
  )

  # Run new function
  res <- lapply(
    seq_len(nrow(grid)),
    function(x) {
      alternative <- grid$alternative[x]
      pValueFUN <- get(grid$pValueFUN[x])
      pValueFUN_args <- list(
        heterogeneity = grid$heterogeneity[x],
        phi = if (grid$heterogeneity[x] == "multiplicative") phi else NULL,
        tau2 = if (grid$heterogeneity[x] == "additive") tau2 else NULL
      )
      hMeanChiSqCI(
        thetahat = thetahat,
        se = se,
        level = level,
        alternative = alternative,
        pValueFUN = pValueFUN,
        pValueFUN_args = pValueFUN_args
      )
    }
  )

  # Get the old function, vectorise it, and run the same inputs
  old_fun <- get_old_FUN(
    path = "https://raw.githubusercontent.com/felix-hof/hMean/main/R/hMeanChiSqCI.R",
    fun_name = "hMeanChiSqCI"
  )

  old_res <- lapply(
    seq_len(nrow(grid)),
    function(x) {
      alternative <- grid$alternative[x]
      pValueFUN <- get(grid$pValueFUN[x])
      pValueFUN_args <- list(
        heterogeneity = grid$heterogeneity[x],
        phi = if (grid$heterogeneity[x] == "multiplicative") phi else NULL,
        tau2 = if (grid$heterogeneity[x] == "additive") tau2 else NULL
      )
      old_fun(
        thetahat = thetahat,
        se = se,
        level = level,
        alternative = alternative,
        pValueFUN = pValueFUN,
        pValueFUN_args = pValueFUN_args
      )
    }
  )

  # compare results
  expect_equal(res, old_res)
})
