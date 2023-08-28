# source globals
source("globals.R")

test_that("Results of pPearsonMu are consistent.", {
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
    stringsAsFactors = FALSE
  )
  # Run new function
  res <- lapply(
    seq_len(nrow(grid)),
    function(x) {
      het <- grid$heterogeneity[x]
      phi <- if (het == "multiplicative") phi else NULL
      tau2 <- if (het == "additive") tau2 else NULL
      pPearsonMu(
        thetahat = thetahat,
        se = se,
        mu = mu,
        phi = phi,
        tau2 = tau2,
        heterogeneity = het,
        alternative = "none",
        check_inputs = TRUE
      )
    }
  )

  # Get the old function, vectorise it, and run the same inputs
  old_fun <- get_old_FUN(
    path = "https://raw.githubusercontent.com/felix-hof/confMeta/main/R/pearson.R",
    fun_name = "pPearsonMu"
  )

  old_vec <- function(
    thetahat,
    se,
    mu,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = c("none", "additive", "multiplicative"),
    alternative = "none",
    check_inputs = TRUE
  ) {
    out <- numeric(length(mu))
    for (i in seq_along(mu)) {
      out[i] <- old_fun(
        thetahat = thetahat,
        se = se,
        mu = mu[i],
        phi = phi,
        tau2 = tau2,
        heterogeneity = heterogeneity,
        alternative = alternative,
        check_inputs = check_inputs
      )
    }
    out
  }
  
  old_res <- lapply(
    seq_len(nrow(grid)),
    function(x) {
      het <- grid$heterogeneity[x]
      phi <- if (het == "multiplicative") phi else NULL
      tau2 <- if (het == "additive") tau2 else NULL
      old_vec(
        thetahat = thetahat,
        se = se,
        mu = mu,
        phi = phi,
        tau2 = tau2,
        heterogeneity = het,
        alternative = "none",
        check_inputs = TRUE
      )
    }
  )
  
  # compare results
  expect_equal(res, old_res)
})
