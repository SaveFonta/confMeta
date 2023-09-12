#' @rdname p_value_functions
#' @order 6
#'
#' @export
#'
#' @examples
#'     # Using Tippett's method to calculate the combined p-value
#'     # for each of the means with multiplicative adjustement for SEs
#'     p_tippett(
#'         estimates = estimates,
#'         SEs = SEs,
#'         mu = mu,
#'         heterogeneity = "multiplicative",
#'         phi = phi
#'     )
p_tippett <- function(
    estimates,
    SEs,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    check_inputs = TRUE
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
    se <- adjust_se(
      se = se,
      heterogeneity = heterogeneity,
      phi = phi,
      tau2 = tau2
    )

    # Get length
    n <- length(estimates)

    # get the z-values
    z <- get_z(estimates = estimates, SEs = SEs, mu = mu)
    # convert them to p-values
    # p <- ReplicationSuccess::z2p(z, "two.sided")
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE) # faster than above
    # Calculate the Tippett statistic
    S_t <- apply(p, 2L, min)
    # Calculate the p-value using the beta distribution
    p <- stats::pbeta(q = S_t, shape1 = 1, shape2 = n, lower.tail = TRUE)
    # return
    return(p)
}
