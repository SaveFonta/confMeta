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
#' The Fisher test statistic for \eqn{n} independent studies is defined as:
#' \deqn{f = -2 \sum_{i=1}^n \log(p_i)}
#' 
#' Under the global null hypothesis, each \eqn{p_i} is assumed to be 
#' uniformly distributed on \eqn{[0, 1]}. The test statistic \eqn{f} therefore follows a 
#' chi-squared distribution with \eqn{2n} degrees of freedom: \eqn{\chi^2_{2n}}. 
#' The combined \emph{p}-value, \eqn{p_F}, is calculated as the probability of observing a 
#' value strictly greater than \eqn{f} from this distribution:
#' \deqn{p_F = \Pr(\chi^2_{2n} > f)}
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
#' @inheritSection p_tippett Best-of-k adjustment
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
#' # Best-of-k adjustment: each study is the most significant of k = 3 experiments
#'  p_fisher(
#'      estimates = estimates,
#'      SEs = SEs,
#'      mu = mu,
#'      heterogeneity = "none",
#'      output_p = "two.sided",
#'      input_p = "greater",
#'      rep = (3,n)
#' )
p_fisher <- function(
    estimates,
    SEs,
    mu = 0,
    heterogeneity = "none",
    phi = NULL,
    tau2 = NULL,
    check_inputs = TRUE,
    input_p = "greater",
    output_p = "two.sided",
    k = rep(1, length(estimates))
) {

  
  # Obtain the matrix of p values of dimension (n_studies x n_mu)
  p <- body_p_value_fun(estimates = estimates,
                        SEs = SEs,
                        mu = mu,
                        heterogeneity = heterogeneity,
                        phi = phi,
                        tau2 = tau2,
                        check_inputs = check_inputs,
                        input_p = input_p,
                        output_p = output_p,
                        k = k)
  
  
  # Get length
  n <- length(estimates) # same as doing  n <- nrow(p)
  
    
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
