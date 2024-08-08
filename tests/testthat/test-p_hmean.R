# source globals
source("globals.R")

test_that("Results of p_hmean are consistent.", {
    # Set inputs
    set.seed(42)
    n <- 15
    estimates <- rnorm(n)
    SEs <- rgamma(n, 5, 5)
    mu <- seq(
        min(estimates) - 0.5 * max(SEs),
        max(estimates) + 0.5 * max(SEs),
        length.out = 100
    )
    phi <- estimate_phi(estimates = estimates, SEs = SEs)
    tau2 <- estimate_tau2(estimates = estimates, SEs = SEs)
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
            p_hmean(
                estimates = estimates,
                SEs = SEs,
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
        path =
        "https://raw.githubusercontent.com/felix-hof/confMeta/main/R/pfun_hmean.R",
        fun_name = "p_hmean"
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
                estimates = estimates,
                SEs = SEs,
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
