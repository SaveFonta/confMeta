#' @title \emph{p}-value functions
#' @rdname p_value_functions
#' @order 1
#'
#' @template estimates
#' @template SEs
#' @template mu
#' @template heterogeneity
#' @template phi
#' @template tau2
#' @template alternative
#' @template check_inputs
#' @param approx Must be either TRUE (default) or FALSE. If TRUE, the p-value is
#'     computed using the normal approximation of the Irwin-Hall distribution
#'     whenever \code{length(estimates) >= 12}. This avoids issues that can lead
#'     to overflow of the double precision floating point numbers R uses for
#'     numeric vectors.
#' @param input_p Either \code{"two.sided"}, \code{"less"}, or \code{"greater"}
#'     (default). Specifies whether two-sided or one-sided p-values are used as
#'     inputs for a p-value combination method.
#'
#' @description
#'     These functions combine individual effect estimates and the corresponding
#'     standard errors into a single \emph{p}-value. Under the hood, all of the
#'     functions transform the estimates and standard errors into \emph{z}, and
#'     subsequently into \emph{p}-values. The resulting \emph{p}-values are
#'     combined into the chosen statistic and an appropriate distribution is
#'     used to derive the combined \emph{p}-value.
#'
#'     All of the \emph{p}-value functions are vectorized over the \code{mu}
#'     argument.
#'
#' @return The corresponding p-values given \eqn{mu} under the null-hypothesis.
#'
#' @note Add references to p-value statistics.
#'
#' @importFrom ReplicationSuccess z2p
#' @importFrom stats pnorm dnorm
#'
#' @export
#'
#' @examples
#'     # Simulating estimates and standard errors
#'     n <- 15
#'     estimates <- rnorm(n)
#'     SEs <- rgamma(n, 5, 5)
#'
#'     # Calculate the between-study variance tau2
#'     tau2 <- estimate_tau2(estimates = estimates, SEs = SEs)
#'     phi <- estimate_phi(estimates = estimates, SEs = SEs)
#'
#'     # Set up a vector of means under the null hypothesis
#'     mu <- seq(
#'       min(estimates) - 0.5 * max(SEs),
#'       max(estimates) + 0.5 * max(SEs),
#'       length.out = 1e5
#'     )
#'
#'     # Using Edgington's method to calculate the combined p-value
#'     # for each of the means with additive adjustement for SEs
#'     p_edgington(
#'         estimates = estimates,
#'         SEs = SEs,
#'         mu = mu,
#'         heterogeneity = "additive",
#'         tau2 = tau2
#'     )
p_edgington <- function(
    estimates,
    SEs,
    mu = 0,
    heterogeneity = "none",
    phi = NULL,
    tau2 = NULL,
    alternative = "two.sided",
    check_inputs = TRUE,
    approx = TRUE,
    input_p = "greater"
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
        check_alternative_arg_edg(alternative = alternative)
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
        ## # p <- ReplicationSuccess::z2p(z, "two.sided")
        p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE) # faster than above
    } else if (input_p == "greater") {
        p <- stats::pnorm(q = z, lower.tail = FALSE)
    } else {
        p <- stats::pnorm(q = z, lower.tail = TRUE)
    }
    # sum up the p-values and calculate the probability
    sp <- pirwinhall(q = colSums(p), n = n, approx = approx)
    if (input_p != "two.sided") {
        sp <- 2*pmin(sp, 1 - sp)
    }
    return(sp)
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
        out <- stats::pnorm(
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
        out[idx_approx] <- stats::pnorm(
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
        out <- stats::dnorm(
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
        out <- vector("double", l_n)
        # fill it with the respective probabilities
        n_norm <- n[idx_approx]
        x_norm <- x[idx_approx]
        out[idx_approx] <- stats::dnorm(
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
    out <- vector("double", l_x)
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
