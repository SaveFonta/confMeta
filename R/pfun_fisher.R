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
    check_inputs = TRUE,
    input_p = "two.sided"
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
    if (input_p == "two.sided") {
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    } else if (input_p == "greater") {
        p <- stats::pnorm(q = z, lower.tail = TRUE)
    } else {
        p <- stats::pnorm(q = z, lower.tail = FALSE)
    }
    # sum up the p-values and calculate the probability
    pfis <- stats::pchisq(
        q = -2 * colSums(log(p)),
        df = 2 * n,
        lower.tail = FALSE
        )
    if (input_p != "two.sided") {
        pfis <- 2*pmin(pfis, 1 - pfis)
    }
    return(pfis)
}
