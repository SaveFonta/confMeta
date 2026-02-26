#' @title Stouffer's method
#' @family p-value combination functions
#'
#' @description
#' Stouffer's method for combining \emph{p}-values across studies
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
#' @param alternative ???
#' @param w Numeric vector of weights. ???
#' 
#' @details
#' Explain how it is computed ???
#'
#' @inherit p_tippett return
#'
#' @importFrom stats pnorm
#' @export
#'
#'
#'
#' @references
#' Stouffer SA, et al. *The American Soldier*. Princeton University Press, 1949.
#'   
#' Held L, Hofmann F, Pawel S. A comparison of combined p-value functions for meta-analysis. *Research Synthesis Methods*, 16:758-785, 2025. 
#'  
#' @examples
#' estimates <- c(0.1, 0.2, 0.3)
#' SEs <- c(0.05, 0.05, 0.1)
#' p_stouffer(estimates, SEs, mu = 0, heterogeneity = "none")

p_stouffer <- function(
    estimates,
    SEs,
    mu = 0,
    heterogeneity = "none",
    phi = NULL,
    tau2 = NULL,
    check_inputs = TRUE,
    alternative = "two.sided",
    w = NULL) {
    # check inputs
    if (check_inputs) {
        check_inputs_p_value(
            estimates = estimates,
            SEs = SEs,
            mu = mu,
            heterogeneity = heterogeneity,
            phi = phi,
            tau2 = tau2
        )
    }

    # recycle `se` if needed
    if (length(SEs) == 1L) SEs <- rep(SEs, length(estimates))

    # adjust se based on heterogeneity model
    SEs <- adjust_se(
        SEs = SEs,
        heterogeneity = heterogeneity,
        phi = phi,
        tau2 = tau2
    )

    ## weights
    if (is.null(w)) {
        ## default set weights to 1/SEs so that Stouffer's method corresponds to
        ## meta-analysis
        w <- 1 / SEs
    } else {
        ## custom weights, recycle if needed
        if (length(w) == 1L) w <- rep(w, length(estimates))
    }

    # compute weighted Stouffer's p-value
    z <- get_z(estimates = estimates, SEs = SEs, mu = mu)
    zs <- colSums(w * z) / sqrt(sum(w^2))
    if (alternative == "two.sided") {
        pstouffer <- 2 * stats::pnorm(q = abs(zs), lower.tail = FALSE)
    } else if (alternative == "less") {
        pstouffer <- stats::pnorm(q = zs, lower.tail = TRUE)
    } else {
        pstouffer <- stats::pnorm(q = zs, lower.tail = FALSE)
    }
    return(pstouffer)
}
