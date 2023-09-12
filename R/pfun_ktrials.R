#' @rdname p_value_functions
#' @order 4

#' @export
#'
#' @examples
#'     # Using the k-trials method to calculate the combined p-value
#'     # for each of the means with multiplicative adjustement for SEs
#'     p_ktrials(
#'         estimates = estimates,
#'         SEs = SEs,
#'         mu = mu,
#'         heterogeneity = "multiplicative",
#'         phi = phi
#'     )
p_ktrials <- function(
    estimates,
    SEs,
    mu,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    alternative = "none",
    check_inputs = TRUE
) {

    if (check_inputs) {
        check_inputs_p_value(
            estimates = estimates,
            SEs = SEs,
            mu = mu,
            heterogeneity = heterogeneity,
            phi = phi,
            tau2 = tau2
        )
        check_alternative_arg_ktr(alternative = alternative)
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

    # Get lengths
    n <- length(estimates)

    if (alternative == "none") {
        z <- get_z(estimates = estimates, SEs = SEs, mu = mu)
        # p <- ReplicationSuccess::z2p(z, "two.sided")
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
        res <- apply(p, 2L, max)^n
    } else {
        stop("Invalid argument 'alternative'.")
    }

    return(res)
}
