#' Calculate the p-value using the Fisher test.
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
#' @template check_inputs
#'
#' @return The corresponding p-values given mu under the null-hypothesis.
#' @examples
#' n <- 15
#' thetahat <- rnorm(n)
#' se <- rgamma(n, 5, 5)
#' mu <- seq(
#'   min(thetahat) - 0.5 * max(se),
#'   max(thetahat) + 0.5 * max(se),
#'   length.out = 1e5
#' )
#' heterogeneity <- "multiplicative"
#' phi <- estimatePhi(thetahat = thetahat, se = se)
#' resP <- pPearsonMu(
#'     thetahat = thetahat,
#'     se = se,
#'     mu = mu,
#'     heterogeneity = heterogeneity,
#'     phi = phi
#' )
#' resH <- hMeanChiSqMu(
#'     thetahat = thetahat,
#'     se = se,
#'     mu = mu,
#'     heterogeneity = heterogeneity,
#'     phi = phi
#' )
#' resTR <- kTRMu(
#'     thetahat = thetahat,
#'     se = se,
#'     mu = mu,
#'     heterogeneity = heterogeneity,
#'     phi = phi
#' )
#' resE <- pEdgingtonMu(
#'     thetahat = thetahat,
#'     se = se,
#'     mu = mu,
#'     heterogeneity = heterogeneity,
#'     phi = phi
#' )
#' resF <- pFisherMu(
#'     thetahat = thetahat,
#'     se = se,
#'     mu = mu,
#'     heterogeneity = heterogeneity,
#'     phi = phi
#' )
#' par(las = 1)
#' matplot(
#'   mu,
#'   cbind(resP, resH, resTR, resE, resF),
#'   type = "l",
#'   lty = 1,
#'   lwd = 3,
#'   ylab = "p-value function"
#' )
#' legend("topleft",
#'   col = 1:5,
#'   lwd = 3,
#'   lty = 1,
#'   c("Pearson", "hMean", "kTR", "Edgington", "Fisher")
#' )
#' title(paste("Phi =", as.character(round(phi, 2))))
#'
#' @export
#'
pFisherMu <- function(
    thetahat,
    se,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
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

    # Get length
    n <- length(thetahat)

    # get the z-values
    z <- get_z(thetahat = thetahat, se = se, mu = mu)
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
