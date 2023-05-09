#' p-Value of Pearson combination test
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
#' resP <- resH <- resTR <- numeric()
#' n <- 15
#' thetahat <- rnorm(n)
#' se <- rgamma(n, 5, 5)
#' mu <- seq(
#'   min(thetahat) - 0.5 * max(se),
#'   max(thetahat) + 0.5 * max(se),
#'   length.out = 1e4
#' )
#' phi <- estimatePhi(thetahat = thetahat, se = se)
#' for(i in 1:length(mu)){
#'   resP[i] <- pPearsonMu(
#'     thetahat=thetahat,
#'     se=se,
#'     mu=mu[i],
#'     heterogeneity="multiplicative",
#'     phi=phi
#'   )
#'   resH[i] <- hMeanChiSqMu(
#'     thetahat=thetahat,
#'     se = se,
#'     mu = mu[i],
#'     heterogeneity = "multiplicative",
#'     phi = phi
#'   )
#'   resTR[i] <- kTRMu(
#'     thetahat = thetahat,
#'     se = se,
#'     mu = mu[i],
#'     heterogeneity = "multiplicative",
#'     phi = phi
#'   )
#' }
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

pPearsonMu <- function(
    thetahat,
    se,
    mu,
    phi = NULL,
    tau2 = NULL,
    alternative = "none",
    heterogeneity = c("none", "additive", "multiplicative"),
    check_inputs = TRUE
) {

    # check inputs
    alternative <- match.arg(alternative)
    heterogeneity <- match.arg(heterogeneity)
    if (check_inputs) {
        stopifnot(
            length(se) == length(thetahat) || length(se) == 1L,
            is.null(phi) || is.numeric(phi) && is.finite(phi),
            is.null(tau2) || is.numeric(tau2) && is.finite(tau2),
            length(mu) == 1L
        )
    }

    # recycle `se` if needed
    if (length(se) == 1L) {
        se <- rep(se, length(thetahat))
    }

    # construct denominator
    se <- switch(
        heterogeneity,
        "none" = se,
        "additive" = sqrt(se^2 + tau2),
        "multiplicative" = sqrt(se^2 * phi)
    )

    # Get lengths
    n <- length(thetahat)

    if (alternative == "none") {
        z <- (thetahat - mu) / se
        p <- ReplicationSuccess::z2p(z, alternative = "two.sided")
        tp <- -2 * sum(log(1 - p))
        p <- stats::pchisq(q = tp, df = 2 * n, lower.tail = TRUE)
    }
    return(p)
}
