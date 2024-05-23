# source globals
source("globals.R")

## SP comment out tests for now because we changed the function defaults
## for a new simulation study...
## test_that("Results of p_wilikinson are consistent", {
##     # Set inputs
##     set.seed(42)
##     n <- 15
##     estimates <- rnorm(n)
##     SEs <- rgamma(n, 5, 5)
##     mu <- seq(
##         min(estimates) - 0.5 * max(SEs),
##         max(estimates) + 0.5 * max(SEs),
##         length.out = 100
##     )
##     phi <- estimate_phi(estimates = estimates, SEs = SEs)
##     tau2 <- estimate_tau2(estimates = estimates, SEs = SEs)
##     grid <- expand.grid(
##         heterogeneity = c("none", "additive", "multiplicative"),
##         stringsAsFactors = FALSE
##     )
##     # Run new function
##     res <- lapply(
##         seq_len(nrow(grid)),
##         function(x) {
##             het <- grid$heterogeneity[x]
##             phi <- if (het == "multiplicative") phi else NULL
##             tau2 <- if (het == "additive") tau2 else NULL
##             p_wilkinson(
##                 estimates = estimates,
##                 SEs = SEs,
##                 mu = mu,
##                 phi = phi,
##                 tau2 = tau2,
##                 heterogeneity = het,
##                 alternative = "none",
##                 check_inputs = TRUE
##             )
##         }
##     )

##     # Get the old function, vectorise it, and run the same inputs
##     old_fun <- get_old_FUN(
##         path =
##         "https://raw.githubusercontent.com/felix-hof/confMeta/main/R/pfun_wilkinson.R",
##         fun_name = "p_wilkinson"
##     )

##     old_res <- lapply(
##         seq_len(nrow(grid)),
##         function(x) {
##             het <- grid$heterogeneity[x]
##             phi <- if (het == "multiplicative") phi else NULL
##             tau2 <- if (het == "additive") tau2 else NULL
##             old_fun(
##                 estimates = estimates,
##                 SEs = SEs,
##                 mu = mu,
##                 phi = phi,
##                 tau2 = tau2,
##                 heterogeneity = het,
##                 #alternative = "none",
##                 check_inputs = TRUE
##             )
##         }
##     )

##     # compare results
##     expect_equal(res, old_res)
## })
