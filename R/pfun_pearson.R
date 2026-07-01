#' @title Pearson's method
#' @family p-value combination functions
#'
#' @description
#' Pearson's method for combining \emph{p}-values across studies. 
#' The method transforms the complement of the individual study \emph{p}-values 
#' using the natural logarithm and evaluates the sum against a chi-squared 
#' distribution.
#'
#' @inheritParams p_tippett
#'
#' @details
#' The Pearson test statistic for \eqn{n} independent studies is defined as:
#' \deqn{g = -2 \sum_{i=1}^n \log(1 - p_i)}
#' 
#' Under the global null hypothesis, each \eqn{p_i} is assumed to be 
#' uniformly distributed on \eqn{[0, 1]}. The test statistic \eqn{g} therefore follows a 
#' chi-squared distribution with \eqn{2n} degrees of freedom: \eqn{\chi^2_{2n}}. 
#' The combined \emph{p}-value, \eqn{p_P}, is calculated as the probability of observing a 
#' value less than or equal to \eqn{g} from this distribution:
#' \deqn{p_P = \Pr(\chi^2_{2n} \leq g)}
#'
#' \strong{Important note on orientation:} Unlike Edgington's method, Pearson's method 
#' is \emph{not} orientation-invariant. The combined \emph{p}-value depends on 
#' the direction of the one-sided \emph{p}-values (controlled by the 
#' \code{input_p} argument). 
#' 
#' Specifically, Pearson's and Fisher's methods are mirrored. Computing 
#' the Pearson combined \emph{p}-value for the "greater" alternative is equal to 
#' 1 minus the Fisher combined \emph{p}-value for the "less" alternative.
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
#' Pearson K. On a method of determining whether a sample of size n supposed to 
#' have been drawn from a parent population having a known probability integral 
#' has probably been drawn at random. *Biometrika*, 25:379-410, 1933.
#' \doi{10.2307/2332290}
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
#' # Using Pearson's method to calculate the combined p-value
#' p_pearson(
#'      estimates = estimates,
#'      SEs = SEs,
#'      mu = mu,
#'      heterogeneity = "none",
#'      output_p = "two.sided",
#'      input_p = "greater"
#' )
#' 
#' # Best-of-k adjustment: each study is the most significant of k = 3 experiments
#' p_pearson(
#'      estimates = estimates,
#'      SEs = SEs,
#'      mu = mu,
#'      heterogeneity = "none",
#'      output_p = "two.sided",
#'      input_p = "greater",
#'      k = rep(3,n)
#' )
p_pearson <- function(
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
  
    
    # Pearson exact calculation: P(ChiSq_2k <= -2 * sum(log(1-p)))
    ppearson <- stats::pchisq(
        q = -2 * colSums(log1p(-p)),
        df = 2 * n,
        lower.tail = TRUE
    )
    if (output_p == "two.sided" && input_p != "two.sided") {
      ppearson <- 2 * pmin(ppearson, 1 - ppearson)
    }
    
    return(ppearson)
}


