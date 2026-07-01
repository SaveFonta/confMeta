#' @title Wilkinson's method
#' @family p-value combination functions
#' 
#' @description
#' Wilkinson's method for combining \emph{p}-values across studies. 
#' This method evaluates the maximum individual study \emph{p}-value.
#'
#' @inheritParams p_tippett
#'
#' @details
#' The Wilkinson combined \emph{p}-value for \eqn{n} independent studies is defined as:
#' \deqn{p_W = \max\{p_1, \dots, p_n\}^n} 
#'
#' Under the global null hypothesis, each \eqn{p_i} is assumed to be 
#' uniformly distributed on \eqn{[0, 1]}.
#' 
#' \strong{Important note on orientation:} Unlike Edgington's method, Wilkinson's method 
#' is \emph{not} orientation-invariant. The combined \emph{p}-value depends on 
#' the direction of the one-sided \emph{p}-values (controlled by the 
#' \code{input_p} argument). 
#' 
#' Specifically, Wilkinson's and Tippett's methods are mirrored. Computing 
#' the Wilkinson combined \emph{p}-value for the "greater" alternative is equal to 
#' 1 minus the Tippett combined \emph{p}-value for the "less" alternative.
#'
#' @inheritSection p_tippett Output p-value
#' 
#' @inheritSection p_tippett Best-of-k adjustment
#'
#' @inherit p_tippett return
#'
#' @importFrom stats pnorm
#'
#' @export
#'
#' @references
#' Wilkinson B. A statistical consideration in psychological research. *Psychological Bulletin*, 48(2):156-158, 1951. 
#' \doi{10.1037/h0059111}
#' 
#' Held, L, Hofmann, F, Pawel, S. (2025). A comparison of combined *p*-value
#' functions for meta-analysis. *Research Synthesis Methods*, 16:758-785.
#' \doi{10.1017/rsm.2025.26}
#'
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
#' # Using Wilkinson's method to calculate the combined p-value
#' p_wilkinson(
#'     estimates = estimates,
#'     SEs = SEs,
#'     mu = mu,
#'     heterogeneity = "none",
#'     output_p = "two.sided",
#'     input_p = "greater"
#' )
#' 
#' # Best-of-k adjustment: each study is the most significant of k = 3 experiments
#' p_wilkinson(
#'     estimates = estimates,
#'     SEs = SEs,
#'     mu = mu,
#'     heterogeneity = "none",
#'     output_p = "two.sided",
#'     input_p = "greater",
#'     k = rep(3,n)
#' )
p_wilkinson <- function(
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
  
  
  #max(p)^n
  sp <- apply(p, 2L, max)^n
  
  # Symmetrize to two-sided ONLY if requested and if inputs weren't already two-sided
  if (output_p == "two.sided" && input_p != "two.sided") {
    sp <- 2 * pmin(sp, 1 - sp)
  }
  
  return(sp)
}
