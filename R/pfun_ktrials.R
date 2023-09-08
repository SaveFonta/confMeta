#' Calculate the p-value using the k-Trials method.

#' @details
#' The function is vectorized over the argument \code{mu}.
#'
#' @template thetahat
#' @template se
#' @template mu
#' @template phi
#' @template tau2
#' @template heterogeneity
#' @template alternative
#' @template check_inputs
#'
#' @return The corresponding p-value given mu under the null-hypothesis.
#' @export
#'
#' @examples
#' thetahat <- c(0.38041824, -0.22170681, -0.09155859)
#' se <- c(0.2256202, 0.2432796, 0.1051743)
#' phi <- estimatePhi(thetahat = thetahat, se = se)
#' tau2 <- estimateTau2(
#'     thetahat,
#'     se,
#'     control = list(stepadj = 0.5, maxiter = 1000, threshold = 1e-6)
#' )
#' mymu <- seq(
#'     min(thetahat - 3 * se),
#'     max(thetahat + 3 * se),
#'     length.out = 1000
#' )
#' p_add <- kTRMu(
#'     thetahat, se,
#'     mymu,
#'     tau2 = tau2,
#'     phi = phi,
#'     heterogeneity = "additive"
#' )
#' p_mult <- kTRMu(
#'     thetahat,
#'     se,
#'     mymu,
#'     tau2 = tau2,
#'     phi = phi,
#'     heterogeneity = "multiplicative"
#' )
p_ktrials <- function(
    thetahat,
    se,
    mu,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    alternative = "none",
    check_inputs = TRUE
) {

    if (check_inputs) {
        check_inputs_p_value(
            thetahat = thetahat,
            se = se,
            mu = mu,
            heterogeneity = heterogeneity,
            phi = phi,
            tau2 = tau2
        )
        check_alternative_arg_ktr(alternative = alternative)
    }

    # recycle `se` if needed
    if (length(se) == 1L) se <- rep(se, length(thetahat))

    # adjust se based on heterogeneity model
    se <- adjust_se(
        se = se,
        heterogeneity = heterogeneity,
        phi = phi,
        tau2 = tau2
    )

    # Get lengths
    n <- length(thetahat)

    if (alternative == "none") {
        z <- get_z(thetahat = thetahat, se = se, mu = mu)
        # p <- ReplicationSuccess::z2p(z, "two.sided")
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
        res <- apply(p, 2L, max)^n
    } else {
        stop("Invalid argument 'alternative'.")
    }

    return(res)
}
