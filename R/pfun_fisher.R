#' @rdname p_value_functions
#' @order 2
#'
#' @examples
#'     # Using Fisher's method to calculate the combined \emph{p}-value
#'     # for each of the means with multiplicative adjustement for SEs
#'     p_fisher(
#'         estimates = estimates,
#'         SEs = SEs,
#'         mu = mu,
#'         heterogeneity = "multiplicative",
#'         phi = phi
#'     )
#'
#' @export
#'
p_fisher <- function(
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
    SEs <- adjust_se(
      SEs = SEs,
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
    # sum up the p-values and calculate the probability
    p <- stats::pchisq(
        q = -2 * colSums(log(p)),
        df = 2 * n,
        lower.tail = FALSE
    )
    return(p)
}
