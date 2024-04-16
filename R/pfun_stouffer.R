#' @rdname p_value_functions
#' @order 7
#'
#' @export
#'
#' @examples
#'     # Using weighted Stouffer's method to calculate the combined p-value for
#'     # each of the means with multiplicative adjustement for SEs
#'     p_stouffer(
#'         estimates = estimates,
#'         SEs = SEs,
#'         mu = mu,
#'         heterogeneity = "multiplicative",
#'         phi = phi
#'     )
p_stouffer <- function(
    estimates,
    SEs,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    alternative = "two.sided",
    check_inputs = TRUE,
    w = NULL
) {

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
        w <- 1/SEs
    } else {
        ## custom weights, recycle if needed
        if (length(w) == 1L) w <- rep(w, length(estimates))
    }

    # Get lengths
    n <- length(estimates)

    # compute weighted Stouffer's p-value
    z <- get_z(estimates = estimates, SEs = SEs, mu = mu)
    zs <- colSums(w*z)/sqrt(sum(w^2))
    if (alternative == "two.sided") {
        pstouffer <- 2*stats::pnorm(q = abs(zs), lower.tail = FALSE)
    } else if (alternative == "less") {
        pstouffer <- stats::pnorm(q = zs, lower.tail = TRUE)
    } else {
        pstouffer <- stats::pnorm(q = zs, lower.tail = FALSE)
    }
    return(pstouffer)
}
