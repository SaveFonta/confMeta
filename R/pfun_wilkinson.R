#' @rdname p_value_functions
#' @order 4

#' @export
#'
#' @examples
#'     # Using Wilkinson's method to calculate the combined p-value
#'     # for each of the means with multiplicative adjustement for SEs
#'     p_wilkinson(
#'         estimates = estimates,
#'         SEs = SEs,
#'         mu = mu,
#'         heterogeneity = "multiplicative",
#'         phi = phi
#'     )
p_wilkinson <- function(
    estimates,
    SEs,
    mu,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    alternative = "none",
    check_inputs = TRUE,
    input_p = "greater"
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
        if (input_p == "two.sided") {
            ## # p <- ReplicationSuccess::z2p(z, "two.sided")
            p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE) # faster than above
        } else if (input_p == "greater") {
            p <- stats::pnorm(q = z, lower.tail = FALSE)
        } else {
            p <- stats::pnorm(q = z, lower.tail = TRUE)
        }
        res <- apply(p, 2L, max)^n
        if (input_p != "two.sided") {
            res <- 2*pmin(res, 1 - res)
        }
    } else {
        stop("Invalid argument 'alternative'.")
    }

    return(res)
}
