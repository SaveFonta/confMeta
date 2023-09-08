#' Calculate the p-value using the Pearson combination test.
#'
#' @details
#' The function is is vectorized over the \code{mu} argument.
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
#' @return The corresponding p-values given mu under the null-hypothesis.
#' @export
#'
#' @examples
#' n <- 15
#' thetahat <- rnorm(n)
#' se <- rgamma(n, 5, 5)
#' mu <- seq(
#'   min(thetahat) - 0.5 * max(se),
#'   max(thetahat) + 0.5 * max(se),
#'   length.out = 1e5
#' )
#' phi <- estimatePhi(thetahat = thetahat, se = se)
#' resP <- pPearsonMu(
#'     thetahat = thetahat,
#'     se = se,
#'     mu = mu,
#'     heterogeneity = "multiplicative",
#'     phi = phi
#' )
#' resH <- hMeanChiSqMu(
#'     thetahat = thetahat,
#'     se = se,
#'     mu = mu,
#'     heterogeneity = "multiplicative",
#'     phi = phi
#' )
#' resTR <- kTRMu(
#'     thetahat = thetahat,
#'     se = se,
#'     mu = mu,
#'     heterogeneity = "multiplicative",
#'     phi = phi
#' )
#' par(las=1)
#' matplot(
#'   mu,
#'   cbind(resP, resH, resTR),
#'   type = "l",
#'   lty = 1,
#'   lwd = 2,
#'   ylab = "p-value function"
#' )
#' legend("topleft", col=c(1,2,3), lwd=2, lty=1, c("Pearson", "hMean","kTR"))
#' title(paste("Phi =", as.character(round(phi, 2))))

p_pearson <- function(
    thetahat,
    se,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = c("none", "additive", "multiplicative"),
    alternative = "none",
    check_inputs = TRUE
) {

    # check inputs
    if (check_inputs) {
        check_inputs_p_value(
            thetahat = thetahat,
            se = se,
            mu = mu,
            heterogeneity = heterogeneity,
            phi = phi,
            tau2 = tau2
        )
        check_alternative_arg_pearson(alternative = alternative)
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

    # implement alternatives
    if (alternative == "none") {
        z <- get_z(thetahat = thetahat, se = se, mu = mu)
        # ReplicationSuccess::z2p
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
        tp <- apply(p, 2L, function(x) -2 * sum(log(1 - x)))
        p <- stats::pchisq(q = tp, df = 2 * n, lower.tail = TRUE)
    } else {
        stop("Invalid argument 'alternative'.")
    }
    return(p)
}
