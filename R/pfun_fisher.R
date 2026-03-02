#' @title Fisher's method
#' @family p-value combination functions
#'
#' @description
#' Fisher's method for combining \emph{p}-values across studies. 
#' The method transforms the individual study \emph{p}-values using the natural 
#' logarithm and evaluates the sum against a chi-squared distribution.
#'
#' @inheritParams p_tippett
#'
#' @details
#' The Fisher test statistic for \eqn{k} studies is defined as:
#' \deqn{f = -2 \sum_{i=1}^k \log(p_i)}
#' 
#' Under the global null hypothesis, each \eqn{p_i} is assumed to be 
#' uniformly distributed on \eqn{[0, 1]}. The test statistic \eqn{f} therefore follows a 
#' chi-squared distribution with \eqn{2k} degrees of freedom: \eqn{\chi^2_{2k}}. 
#' The combined \emph{p}-value, \eqn{p_F}, is calculated as the probability of observing a 
#' value strictly greater than \eqn{f} from this distribution:
#' \deqn{p_F = \Pr(\chi^2_{2k} > f)}
#'
#' \strong{Important note on orientation:} Unlike Edgington's method, Fisher's method 
#' is \emph{not} orientation-invariant. The combined \emph{p}-value depends on 
#' the direction of the one-sided \emph{p}-values (controlled by the 
#' \code{input_p} argument). 
#' 
#' Specifically, Fisher's and Pearson's methods are mirrored. Computing 
#' the Fisher combined \emph{p}-value for the "greater" alternative is equal to
#' 1 minus the Pearson combined \emph{p}-value for the "less" alternative.
#'
#' @inheritSection p_tippett Output p-value
#'
#' @inherit p_tippett return
#'
#' @importFrom stats pnorm pchisq
#'
#' @export
#'
#' @references
#' Fisher R.A. *Statistical Methods for Research Workers*. 4th ed. Oliver & Boyd; 1932. 
#' \doi{10.1093/oso/9780198522294.002.0003}
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
#'    min(estimates) - 0.5 * max(SEs),
#'    max(estimates) + 0.5 * max(SEs),
#'    length.out = 100
#' )
#'
#' # Using Fisher's method to calculate the combined p-value
#' p_fisher(
#'      estimates = estimates,
#'      SEs = SEs,
#'      mu = mu,
#'      heterogeneity = "none",
#'      output_p = "two.sided",
#'      input_p = "greater"
#' )
#' 
p_fisher <- function(
    estimates,
    SEs,
    mu = 0,
    heterogeneity = "none",
    phi = NULL,
    tau2 = NULL,
    check_inputs = TRUE,
    input_p = "greater",
    output_p = "two.sided"
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
    
    ## Fisher calculation: P(ChiSq_2k > -2 * sum(log(p)))
    pfis <- stats::pchisq(
        q = -2 * colSums(log(p)),
        df = 2 * n,
        lower.tail = FALSE
    )
    
    if (output_p == "two.sided" && input_p != "two.sided") {
      pfis <- 2 * pmin(pfis, 1 - pfis)
    }
    
    return(pfis)
}
