#' @title Edgington's method
#' @family p-value combination functions
#'
#' @description
#' Edgington’s method for combining \emph{p}-values across studies. 
#' The method forms a sum of individual study \emph{p}-values and evaluates 
#' it against the exact or approximate null distribution.
#'
#' Under the global null hypothesis, the null distribution of the sum is given 
#' by the Irwin–Hall distribution. For a weighted generalization of this 
#' procedure, see \code{\link{p_edgington_w}}. 
#'
#' @inheritParams p_tippett
#' @param approx Logical (default \code{TRUE}). If \code{TRUE}, use a normal
#'     approximation for the sum of \emph{p}-values when \eqn{k\geq 12} to 
#'     avoid numerical overflow issues.
#'
#' @details
#' The classical Edgington statistic is defined for \eqn{k} studies as
#' \deqn{S = \sum_{i=1}^k p_i,}
#' where \eqn{p_i} are individual study \emph{p}-values. Under the global null 
#' hypothesis, each \eqn{p_i} is assumed to be 
#' uniformly distributed on \eqn{[0, 1]}. 
#' 
#' \strong{Important note on orientation:} Edgington's method is orientation-invariant. 
#' The combined \emph{p}-value is symmetric with respect to the direction of the 
#' one-sided \emph{p}-values (controlled by the \code{input_p} argument). 
#' 
#' Specifically, computing the Edgington combined \emph{p}-value for the "greater" 
#' alternative results in 1 minus the Edgington combined \emph{p}-value for 
#' the "less" alternative.
#' 
#' @section Null Distribution and Approximation:
#' The combined \emph{p}-value, \eqn{p_E}, is the probability of observing a sum 
#' less than or equal to \eqn{S} under the null hypothesis. This is computed in 
#' one of two ways:
#' \itemize{
#'   \item \strong{Exact Method:} The function uses the exact Irwin-Hall distribution 
#'     to compute the combined \emph{p}-value:
#'     \deqn{p_E = \frac{1}{k!} \sum_{j=0}^{\lfloor S \rfloor} (-1)^j \binom{k}{j} (S - j)^k}
#'   \item \strong{Normal Approximation:} For a large number of studies (\eqn{k \geq 12}), 
#'     the distribution of the sum is approximated by a Normal distribution with:
#'     \deqn{\mathrm{E}[S] = \frac{k}{2}}
#'     \deqn{\mathrm{Var}(S) = \frac{k}{12}}
#' }
#'
#' @inheritSection p_tippett Output p-value
#'
#' @inherit p_tippett return
#'
#' @importFrom stats pnorm dnorm
#'
#' @export
#'
#' @references
#' Edgington, E. S. (1972). An additive method for combining probability values from
#'   independent experiments. *The Journal of Psychology*, 80(2): 351-363.
#'   \doi{10.1080/00223980.1972.9924813}
#'   
#' Held, L, Hofmann, F, Pawel, S. (2025). A comparison of combined *p*-value
#' functions for meta-analysis. *Research Synthesis Methods*, 16:758-785.
#' \doi{10.1017/rsm.2025.26}
#'
#' @examples
#' # Simulating estimates and standard errors
#' n <- 15
#' estimates <- rnorm(n)
#' SEs <- rgamma(n, 5, 5)
#'
#' # Set up a vector of means under the null hypothesis
#' mu <- seq(
#'   min(estimates) - 0.5 * max(SEs),
#'   max(estimates) + 0.5 * max(SEs),
#'   length.out = 100
#' )
#'
#' # Using Edgington's method to calculate the combined p-value
#' p_edgington(
#'     estimates = estimates,
#'     SEs = SEs,
#'     mu = mu,
#'     heterogeneity = "none",
#'     output_p = "two.sided",
#'     input_p = "greater",
#'     approx = TRUE
#' )
p_edgington <- function(
    estimates,
    SEs,
    mu = 0,
    heterogeneity = "none",
    phi = NULL,
    tau2 = NULL,
    check_inputs = TRUE,
    input_p = "greater",
    output_p = "two.sided",
    approx = TRUE
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
      check_output_p_arg(output_p = output_p)
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
    p <- switch(input_p,
                "two.sided" = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
                "greater"   = stats::pnorm(z, lower.tail = FALSE),
                "less"      = stats::pnorm(z, lower.tail = TRUE),
                stop("input_p must be 'greater','less','two.sided'")
    )
    
    p <- as.matrix(p)
    
    # sum up the p-values and calculate the probability
    sp <- pirwinhall(q = colSums(p), n = n, approx = approx)
    
    # Symmetrize to two-sided ONLY if requested and if inputs weren't already two-sided
    if (output_p == "two.sided" && input_p != "two.sided") {
      sp <- 2 * pmin(sp, 1 - sp)
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

## ??? --> maybe delete this overflow example?

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
        if (l_q == 1L) q <- rep(q, l_n) 
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
#pirwinhall(q = 2, n = c(3, 5, 50), approx = TRUE)

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
        s_l_n <- seq_along(n) #AHA! previous version used to have seq_len(n), but need to use seq_along(n)!! Or seq_len(l_n)
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
#
#pirwinhall_approx(q = c(3, 6.0, 25), n = c(6, 12, 50))

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
#here we can text two SEPARATE cases at once
# Case 1: q = 1.5, n = 3
# Case 2: q = 2.0, n = 4
#pirwinhall_vec(q = c(1.5, 2.0), n = c(3, 4))


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

#Example: we want prob that the sum of three Uniforms is <= 1.5
#pirwinhall1(q = 1.5, n = 3)


###################################################################


#  ??? --> the functions above here are completely useless! I was thinking to remove them completely even though are nice functions
# But they are not used at all, since we only need the CDF (pirwinhall)


###################################################################


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
