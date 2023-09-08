#' @rdname p_value_functions
#' @order 3
#'
#' @template w
#' @template distr
#'
#' @examples
#'     # Using the harmonic mean method to calculate the combined p-value
#'     # for each of the means with additive adjustment for SEs.
#'     p_hmean(
#'         estimates = estimates,
#'         SEs = SEs,
#'         mu = mu,
#'         heterogeneity = "additive",
#'         tau2 = tau2,
#'         distr = "chisq"
#'     )
#'
#' @export
p_hmean <- function(
    thetahat,
    se,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    alternative = "none",
    check_inputs = TRUE,
    w = rep(1, length(thetahat)),
    distr = "chisq"
) {

    # Check inputs
    if (check_inputs) {
        check_inputs_p_value(
            thetahat = thetahat,
            se = se,
            mu = mu,
            heterogeneity = heterogeneity,
            phi = phi,
            tau2 = tau2
        )
        check_alternative_arg_hmean(alternative = alternative)
        check_distr_arg(distr = distr)
        check_w_arg(w = w, thetahat = thetahat)
    }

    # match arguments
    if (length(se) == 1L) se <- rep(se, length(thetahat))

    # adjust se based on heterogeneity model
    se <- adjust_se(
        se = se,
        heterogeneity = heterogeneity,
        phi = phi,
        tau2 = tau2
    )

    # store lengths of input vector
    n <- length(thetahat)

    # Calculate harmonic mean test statistic
    sw <- sum(sqrt(w))^2
    z <- get_z(thetahat = thetahat, se = se, mu = mu)
    zh2 <- apply(z, 2L, function(z) sw / sum(w / z^2))
    # Calculate the p-value
    res <- switch(
        distr,
        "chisq" = stats::pchisq(zh2, df = 1L, lower.tail = FALSE),
        "f" = stats::pf(zh2, df1 = 1L, df2 = n - 1, lower.tail = FALSE)
    )

    if (alternative != "none") {
        check_g <- apply(z, 2L, function(z) min(z) >= 0)
        check_l <- apply(z, 2L, function(z) max(z) <= 0)
        cond <- check_g | check_l
        res <- switch(
            alternative,
            "greater" = ifelse(cond, res / 2^n, NaN),
            "less" = ifelse(cond, res / 2^n, NaN),
            "two.sided" = ifelse(cond, res / 2^(n - 1), NaN),
        )
    }

    # return
    res
}
