#' @title Harmonic mean's method
#' @family p-value combination functions
#' 
#' @description
#' Combines study-level results using the harmonic mean of squared \eqn{z}-statistics.#' 
#' 
#' @param estimates Numeric vector of study-level effect estimates.
#' @param SEs Numeric vector of corresponding standard errors.
#' @param mu Numeric **scalar or vector** of null values for the overall effect
#'      (default: 0). 
#' @param heterogeneity Character string: \code{"none"} (default),
#'      \code{"additive"}, or \code{"multiplicative"}. Determines whether
#'      standard errors are adjusted for between-study heterogeneity using
#'      \code{tau2} or \code{phi}.
#' @param phi Multiplicative heterogeneity parameter (if applicable).
#' @param tau2 Additive heterogeneity parameter (if applicable).
#' @param check_inputs Logical (default \code{TRUE}). If \code{TRUE},
#'      perform input validation.
#' @param alternative  ???
#' @param w Numeric vector of weights (default: equal weights).
#' @param distr Character string specifying the null distribution: 
#'      \code{"chisq"} (default) or \code{"f"}.
#' 
#' 
#' @details
#' Explain how it is computed ???
#' 
#' 
#' 
#' @inherit p_tippett return
#'
#' @export
#' @importFrom stats pf pchisq
#'
#'
#' @references
#' 
#' ???
#' 
#' Held L, Hofmann F, Pawel S. A comparison of combined p-value functions for meta-analysis. *Research Synthesis Methods*, 16:758-785, 2025.
#' 
#'  
#' @examples
#' estimates <- c(0.5, 0.8, 0.3)
#' SEs <- c(0.1, 0.2, 0.1)
#' p_hmean(estimates, SEs, mu = 0, distr = "f")
#' 

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
