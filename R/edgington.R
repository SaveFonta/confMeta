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
#' heterogeneity <- "none"
#' phi <- NULL
#' # phi <- estimatePhi(thetahat = thetahat, se = se)
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
#' resE <- edgingtonMu(
#'     thetahat = thetahat,
#'     se = se,
#'     mu = mu,
#'     heterogeneity = heterogeneity,
#'     phi = phi
#' )
#' par(las=1)
#' matplot(
#'   mu,
#'   cbind(resP, resH, resTR, resE),
#'   type = "l",
#'   lty = 1,
#'   lwd = 2,
#'   ylab = "p-value function"
#' )
#' legend("topleft",
#'   col=c(1, 2, 3, 4),
#'   lwd=2,
#'   lty=1,
#'   c("Pearson", "hMean","kTR","Edgington")
#' )
#' title(paste("Phi =", as.character(round(phi, 2))))

edgingtonMu <- function(
    thetahat,
    se,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = c("none", "additive", "multiplicative"),
    alternative = c("two.sided", "one.sided"),
    check_inputs = TRUE,
    approx = FALSE
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
    z <- vapply(
        mu,
        function(mu) {
            (thetahat - mu) / se
        },
        double(n)
    )
    if (is.null(dim(z)))
        dim(z) <- c(1L, n)
    # convert them to p-values
    # p <- ReplicationSuccess::z2p(z, "two.sided")
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE) # faster than above
    Fx <- pirwinhall(colSums(p), n = n, approx = approx)
    if (alternative == "two.sided") { # two.sided
        p <- 2 * apply(matrix(c(Fx, 1 - Fx), ncol = 2L), 1L, min)
    } else { # one.sided
        p <- Fx
    }
    return(p)
}

# # Note that at around n = 50, dirwinhall starts to overflow
# n <- 50
# x <- rowSums(matrix(runif(n * 1e5), ncol = n))
# hist(x, prob = TRUE, nclass = 100, ylim = c(0, 0.5))
# x_line <- seq(min(x), max(x), length.out = 1e4)
# lines(x_line, dnorm(x = x_line, mean = n / 2, sd = sqrt(n / 12)), col = 2L)
# lines(x_line, dirwinhall(x = x_line, n = n), col = 3L)
# legend(
#     x = "topright",
#     legend = c("dnorm", "dirwinhall"),
#     col = c(2, 3),
#     lty = 1
# )

################################################################################
# Irwin Hall ===================================================================
################################################################################

## Distribution function of Irwin-Hall distribution ============================

# vectorize over x and n and add some other arguments
pirwinhall <- function(q, n, lower.tail = TRUE, log.p = FALSE, approx = FALSE) {

    # Check arguments
    l_n <- length(n)
    l_q <- length(q)

    repeat_arg <- l_n != l_q
    case_ok <- !repeat_arg || l_n == 1L || l_q == 1L
    if (!case_ok) stop("Invalid argument configuration.")

    # repeat n and x if necessary
    if (repeat_arg) {
        if (l_n == 1L) n <- rep(n, l_q)
        if (l_q == 1L) x <- rep(x, l_n)
    }

    # call pirwinhall
    if (approx) {
        out <- pirwinhall_approx(q, n)
    } else {
        out <- pirwinhall_vec(q, n)
    }

    if (!lower.tail) out <- 1 - out
    if (log.p) log(out) else out
}

# this function uses a normal approximation if n >= 12
# in order to mitigate overflow problems of pirwinhall1
# when n is large
pirwinhall_approx <- function(q, n) {
    # which elements use the normal approximation?
    norm_approx <- n > 11
    n_norm_approx <- sum(norm_approx)
    l_n <- length(n)
    if (n_norm_approx == 0L) {
        # if all(n <= 11), just use pirwinhall_vec for all elements
        out <- pirwinhall_vec(q = q, n = n)
    } else if (n_norm_approx == l_n) {
        # if all(n > 11), just use pnorm for all elements
        out <- pnorm(
            q = q,
            mean = n / 2,
            sd = sqrt(n / 12),
            log.p = FALSE,
            lower.tail = TRUE
        )
    } else {
        # for elements where n <= 11, use pirwinhall_vec, and
        # for all others use pnorm

        # get the indices of elements that use the approximation
        # and of those elements that use regular irwin-hall
        s_l_n <- seq_len(n)
        idx_approx <- s_l_n[norm_approx]
        idx_non_approx <- s_l_n[!norm_approx]
        # create the output vector
        out <- vector("numeric", l_n)
        # fill it with the respective probabilities
        n_norm <- n[idx_approx]
        q_norm <- q[idx_approx]
        out[idx_approx] <- pnorm(
            q = q_norm,
            mean = n_norm / 2,
            sd = sqrt(n_norm / 12),
            log.p = FALSE,
            lower.tail = TRUE
        )
        n_irwhall <- n[idx_non_approx]
        q_irwhall <- q[idx_non_approx]
        out[idx_non_approx] <- pirwinhall_vec(
            q = q_irwhall,
            n = n_irwhall
        )
    }
    # return
    out
}

# vectorized version of pirwinhall1,
# needs x and n to be of equal length
pirwinhall_vec <- function(q, n) {

    l_q <- length(q)
    out <- vector("numeric", l_q)
    for (i in seq_len(l_q)) {
        out[i] <- pirwinhall1(q = q[i], n = n[i])
    }
    out
}

# function for one x (i.e. length(x) == 1)
pirwinhall1 <- function(q, n) {

    if (q <= 0) {
        0
    } else if (q >= n) {
        1
    } else {
        k <- 0:floor(q)
        sum((-1)^k * choose(n = n, k = k) * (q - k)^n) / factorial(n)
    }
}

## Density function of Irwin-Hall distribution =================================

# vectorize over x and n and add some other arguments
dirwinhall <- function(x, n, log = FALSE, approx = FALSE) {

    # Check arguments
    l_n <- length(n)
    l_x <- length(x)

    repeat_arg <- l_n != l_x
    case_ok <- !repeat_arg || l_n == 1L || l_x == 1L
    if (!case_ok) stop("Invalid argument configuration.")

    # repeat n and x if necessary
    if (repeat_arg) {
        if (l_n == 1L) n <- rep(n, l_x)
        if (l_x == 1L) x <- rep(x, l_n)
    }

    # call pirwinhall
    if (approx) {
        out <- dirwinhall_approx(x, n)
    } else {
        out <- dirwinhall_vec(x, n)
    }

    if (log) log(out) else out
}

# this function uses a normal approximation if n >= 12
# in order to mitigate overflow problems of dirwinhall1
# when n is large
dirwinhall_approx <- function(x, n) {
    # which elements use the normal approximation?
    norm_approx <- n > 11
    n_norm_approx <- sum(norm_approx)
    l_n <- length(n)
    if (n_norm_approx == 0L) {
        # if all(n <= 11), just use dirwinhall_vec for all elements
        out <- dirwinhall_vec(x = x, n = n)
    } else if (n_norm_approx == l_n) {
        # if all(n > 11), just use pnorm for all elements
        out <- dnorm(
            x = x,
            mean = n / 2,
            sd = sqrt(n / 12),
            log = FALSE
        )
    } else {
        # for elements where n <= 11, use dirwinhall_vec, and
        # for all others use dnorm

        # get the indices of elements that use the approximation
        # and of those elements that use regular irwin-hall
        s_l_n <- seq_len(n)
        idx_approx <- s_l_n[norm_approx]
        idx_non_approx <- s_l_n[!norm_approx]
        # create the output vector
        out <- vector("numeric", l_n)
        # fill it with the respective probabilities
        n_norm <- n[idx_approx]
        x_norm <- x[idx_approx]
        out[idx_approx] <- dnorm(
            x = x_norm,
            mean = n_norm / 2,
            sd = sqrt(n_norm / 12),
            log = FALSE
        )
        n_irwhall <- n[idx_non_approx]
        x_irwhall <- x[idx_non_approx]
        out[idx_non_approx] <- dirwinhall_vec(
            x = x_irwhall,
            n = n_irwhall
        )
    }
    # return
    out
}

# vectorized version of dirwinhall1,
# needs x and n to be of equal length
dirwinhall_vec <- function(x, n) {

    l_x <- length(x)
    out <- vector("numeric", l_x)
    for (i in seq_len(l_x)) {
        out[i] <- dirwinhall1(x = x[i], n = n[i])
    }
    out
}

# function for one x
dirwinhall1 <- function(x, n) {

    if (x <= 0 || x >= n) {
        0
    } else {
        k <- 0:floor(x)
        sum((-1)^k * choose(n = n, k = k) * (x - k)^(n - 1)) / factorial(n - 1)
    }
}
