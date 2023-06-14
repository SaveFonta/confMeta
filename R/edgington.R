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
# n <- 15
# thetahat <- rnorm(n)
# se <- rgamma(n, 5, 5)
# mu <- seq(
#   min(thetahat) - 0.5 * max(se),
#   max(thetahat) + 0.5 * max(se),
#   length.out = 1e5
# )
# phi <- estimatePhi(thetahat = thetahat, se = se)
# resP <- pPearsonMu(
#     thetahat = thetahat,
#     se = se,
#     mu = mu,
#     heterogeneity = "multiplicative",
#     phi = phi
# )
# resH <- hMeanChiSqMu(
#     thetahat = thetahat,
#     se = se,
#     mu = mu,
#     heterogeneity = "multiplicative",
#     phi = phi
# )
# resTR <- kTRMu(
#     thetahat = thetahat,
#     se = se,
#     mu = mu,
#     heterogeneity = "multiplicative",
#     phi = phi
# )
# resE <- edgingtonMu(
#     thetahat = thetahat,
#     se = se,
#     mu = mu,
#     heterogeneity = heterogeneity,
#     phi = phi
# )
# par(las=1)
# matplot(
#   mu,
#   cbind(resP, resH, resTR),
#   type = "l",
#   lty = 1,
#   lwd = 2,
#   ylab = "p-value function"
# )
# legend("topleft", col=c(1,2,3), lwd=2, lty=1, c("Pearson", "hMean","kTR"))
# title(paste("Phi =", as.character(round(phi, 2))))

edgingtonMu <- function(
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
        alternative = alternative,
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

    # implement alternatives
    if (alternative == "none") {
        z <- vapply(
            mu,
            function(mu) {
                (thetahat - mu) / se
            },
            double(length(thetahat))
        )
        if (is.null(dim(z)))
            dim(z) <- c(1L, length(z))
        # ReplicationSuccess::z2p
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
        fx <- pirwinhall(colSums(p), n = n)
        if (TRUE) { # two.sided
            p <- min(fx, 1 - fx)
        } else { # one.sided
            p <- fx
        }
    } else {
        stop("Invalid argument 'alternative'.")
    }
    return(p)
}



# x <- rowSums(matrix(runif(12000), ncol = 12))
# x <- 28
# n <- 30
# pirwinhall1(x = x, n = n)
# pirwinhall1_log(x = x, n = n)

pirwinhall <- function(x, n, lower.tail = TRUE, log.p = FALSE) {

    l_n <- length(n)
    l_x <- length(x)

    if (l_n != l_x && l_n != 1L)
        stop("Argument `n` must be either length 1 or length(x).")

    if (l_n == 1L && l_x != 1L) n <- rep(n, l_x)

    out <- vapply(
        seq_along(x),
        function(i) {
            pirwinhall1(x = x[i], n = n[i])
        },
        double(1L)
    )

    if (!lower.tail) out <- 1 - out
    if (log.p) log(out) else out
}


# pirwinhall1_log <- function(x, n) {
#
#     if (x <= 0) {
#         0
#     } else if (x >= n) {
#         1
#     } else {
#         denom  <- lfactorial(n)
#         k <- 0:floor(x)
#         log_terms <- vapply(
#             k,
#             function(k) {
#                 lchoose(n, k) + n * log(x - k) - denom
#             },
#             double(1L)
#         )
#         if (0) "null" else if (1) "one"
#         pm <- ifelse((k %% 2), 1, -1)
#         sum(pm * exp(log_terms))
#     }
# }


pirwinhall1 <- function(x, n) {

    if (x <= 0) {
        0
    } else if (x >= n) {
        1
    } else {
        s <- sum(
            vapply(
                0:floor(x),
                function(k) {
                    out <- choose(n = n, k = k) * (x - k)^n
                    if (k %% 2L == 0L) out else -out
                },
                double(1L)
            )
        )
        if (n == 1L) s else s / factorial(n)
    }
}
