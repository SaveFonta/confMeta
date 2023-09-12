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
    estimates,
    SEs,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    alternative = "none",
    check_inputs = TRUE,
    w = rep(1, length(estimates)),
    distr = "chisq"
) {

    # Check inputs
    if (check_inputs) {
        check_inputs_p_value(
            estimates = estimates,
            SEs = SEs,
            mu = mu,
            heterogeneity = heterogeneity,
            phi = phi,
            tau2 = tau2
        )
        check_alternative_arg_hmean(alternative = alternative)
        check_distr_arg(distr = distr)
        check_w_arg(w = w, estimates = estimates)
    }

    # match arguments
    if (length(SEs) == 1L) SEs <- rep(SEs, length(estimates))

    # adjust se based on heterogeneity model
    SEs <- adjust_se(
        SEs = SEs,
        heterogeneity = heterogeneity,
        phi = phi,
        tau2 = tau2
    )

    # store lengths of input vector
    n <- length(estimates)

    # Calculate harmonic mean test statistic
    sw <- sum(sqrt(w))^2
    z <- get_z(estimates = estimates, SEs = SEs, mu = mu)
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
