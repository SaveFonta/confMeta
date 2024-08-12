#' @rdname p_value_functions
#' @order 5
#'
#' @importFrom stats pnorm pchisq
#'
#' @export
#'
#' @examples
#' # Using Pearson's method to calculate the combined p-value
#' # for each of the means with multiplicative adjustement for SEs
#' p_pearson(
#'     estimates = estimates,
#'     SEs = SEs,
#'     mu = mu,
#'     heterogeneity = "multiplicative",
#'     phi = phi
#' )
p_pearson <- function(
    estimates,
    SEs,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    check_inputs = TRUE,
    input_p = "greater") {
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

    # Get lengths
    n <- length(estimates)

    # implement alternatives
    z <- get_z(estimates = estimates, SEs = SEs, mu = mu)
    if (input_p == "two.sided") {
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    } else if (input_p == "greater") {
        p <- stats::pnorm(q = z, lower.tail = FALSE)
    } else {
        p <- stats::pnorm(q = z, lower.tail = TRUE)
    }
    ## # ReplicationSuccess::z2p
    ## p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    ## tp <- apply(p, 2L, function(x) -2 * sum(log(1 - x)))
    ## p <- stats::pchisq(q = tp, df = 2 * n, lower.tail = TRUE)
    ppearson <- stats::pchisq(
        q = -2 * colSums(log(1 - p)),
        df = 2 * n,
        lower.tail = TRUE
    )
    if (input_p != "two.sided") {
        ppearson <- 2 * pmin(ppearson, 1 - ppearson)
    }
    return(ppearson)
}
