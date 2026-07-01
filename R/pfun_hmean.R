#' @title Harmonic mean method
#' @family p-value combination functions
#' 
#' @description
#' Combines study-level results using the harmonic mean combination method
#' 
#' @inheritParams p_tippett
#' @template w
#' @template distr
#' 
#' 
#' @details
#' The harmonic mean statistic for \eqn{n} independent studies is defined as
#' \deqn{x^2 = \frac{(\sum_{i=1}^n \sqrt{w_i})^2}{\sum_{i=1}^n w_i/z_i^2},}
#' where \eqn{z_i} and \eqn{{w_i}} are the individual study \emph{z}-values and weights,
#' respectively. Under the global null hypothesis, each \eqn{z_i} is assumed to
#' follow a standard normal distribution. The harmonic mean statistic then
#' follows a chi-squared distribution with one degree of freedom. The combined
#' \emph{p}-value is the probability of observing a value equal to or
#' greater than \eqn{x^2}: \deqn{p_H = \Pr(\chi^2_{1} > x^2).}
#'
#'
#' @section Best-of-k adjustment:
#' The harmonic mean method does \strong{not} support the best-of-k adjustment.
#' The statistic depends on the \emph{z}-values only through \eqn{z_i^2}, which
#' discards their sign and hence the direction of selection. The \code{k} argument is therefore accepted only for interface
#' consistency with the other combination functions: it must equal 1 for every
#' study (the default).
#' 
#' @inherit p_tippett return
#'
#' @export
#' @importFrom stats pf pchisq
#'
#'
#' @references
#' 
#' Held, L. (2020). The harmonic mean chi-squared test to substantiate
#' scientific findings. *Journal of the Royal Statistical Society: Series C
#' (Applied Statistics)*, 69:697-708. \doi{10.1111/rssc.12410}
#' 
#' Held, L, Hofmann, F, Pawel, S. (2025). A comparison of combined *p*-value
#' functions for meta-analysis. *Research Synthesis Methods*, 16:758-785.
#' \doi{10.1017/rsm.2025.26}
#' 
#'  
#' @examples
#' estimates <- c(0.5, 0.8, 0.3)
#' SEs <- c(0.1, 0.2, 0.1)
#' p_hmean(estimates, SEs, mu = 0, distr = "f")
#' 
p_hmean <- function(
    estimates,
    SEs,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    check_inputs = TRUE,
    w = rep(1, length(estimates)),
    distr = "chisq",
    k = rep(1, length(estimates)) #useless, just for compatibility with make_function in confMeta
) {

    # Check inputs
    if (check_inputs) {
        check_inputs_p_value(
            estimates = estimates,
            SEs = SEs,
            mu = mu,
            heterogeneity = heterogeneity,
            phi = phi,
            tau2 = tau2
        )
        check_distr_arg(distr = distr)
        check_w_arg(w = w, estimates = estimates)
        check_k_arg(k = k, estimates = estimates) 
    }

    # match arguments
    if (length(SEs) == 1L) SEs <- rep(SEs, length(estimates))

    # adjust se based on heterogeneity model
    SEs <- adjust_se(
        SEs = SEs,
        heterogeneity = heterogeneity,
        phi = phi,
        tau2 = tau2
    )

    # store lengths of input vector
    n <- length(estimates)

    # Calculate harmonic mean test statistic
    sw <- sum(sqrt(w))^2
    z <- get_z(estimates = estimates, SEs = SEs, mu = mu)
    
    
    # I am not sure how to specify the k, since here the Z scores are squared and we lose 
    # directionality.... for the moment comment out and stop in case using k different than 1
    if (any(k > 1L)) {
      stop(
        "Harmonic Mean method is non-directional (z^2) and doesn't support the best-of-k adjustment. ",
        call. = FALSE
      )
    }
    
    ## if (!all(k == 1L)) {
    ##   z <- adjust_z_for_best_k(z = z, k = k, input_p = "greater")
    ## }
    
    
    zh2 <- sw / colSums(w / z^2)
   
    # Calculate the p-value
    res <- switch(
        distr,
        "chisq" = stats::pchisq(zh2, df = 1L, lower.tail = FALSE),
        "f" = stats::pf(zh2, df1 = 1L, df2 = n - 1, lower.tail = FALSE)
    )

    # return
    res
}
