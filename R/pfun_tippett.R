#' @title Tippett's method
#' @family p-value combination functions
#'
#' @description
#' Tippett's method for combining \emph{p}-values across studies. 
#' This method evaluates the minimum individual study \emph{p}-value.
#'
#' @template estimates
#' @template SEs
#' @template mu 
#' @template heterogeneity 
#' @template phi 
#' @template tau2
#' @template check_inputs
#' @template output_p
#' @param input_p Type of study-level \emph{p}-values used in the combination:
#'     \code{"greater"} (default), \code{"less"}, or \code{"two.sided"}.
#'     If \code{"greater"} or \code{"less"}, one-sided \emph{p}-values are
#'     combined.
#' @template k     
#'
#' @details
#' The Tippett combined \emph{p}-value, \eqn{p_T}, for \eqn{n} independent studies is 
#' defined directly as:
#' \deqn{p_T = 1 - (1 - \min\{p_1, \dots, p_k\})^n}
#'
#' Under the global null hypothesis, each \eqn{p_i} is assumed to be 
#' uniformly distributed on \eqn{[0, 1]}.
#' 
#' \strong{Important note on orientation:} Unlike Edgington's method, Tippett's method 
#' is \emph{not} orientation-invariant. The combined \emph{p}-value depends on 
#' the direction of the one-sided \emph{p}-values (controlled by the 
#' \code{input_p} argument). 
#' 
#' Specifically, Tippett's and Wilkinson's methods are mirrored. Computing 
#' the Tippett combined \emph{p}-value for the "greater" alternative is equal to
#' 1 minus the Wilkinson combined \emph{p}-value for the "less" alternative.
#'
#' @section Output p-value:
#' The final output depends on the \code{output_p} and \code{input_p} arguments:
#' \itemize{
#'   \item If \code{output_p = "two.sided"} (the default) and the inputs are 
#'     one-sided (\code{input_p} is \code{"greater"} or \code{"less"}), the function 
#'     combines the one-sided \emph{p}-values to obtain the intermediate combined 
#'     \emph{p}-value \eqn{p_c}, and returns a \strong{symmetrized, two-sided \emph{p}-value}:
#'     \eqn{p_{2s} = 2 \min(p_c, 1 - p_c)}.
#'   \item If \code{output_p = "one.sided"}, the function 
#'     returns the inherently one-sided combined \emph{p}-value \eqn{p_c} 
#'     \strong{directly, without symmetrization}.
#'   \item If \code{input_p} is \code{"two.sided"}, the input \eqn{p_i}
#'     are already two-sided, and no further symmetrization is applied. 
#' }
#' 
#' @section Best-of-k adjustment:
#' The \code{k} argument adjusts for selective reporting, where the
#' reported result of study \eqn{i} is the most significant out of
#' \eqn{k_i} independent experiments. Under the global null hypothesis the
#' smallest of \eqn{k_i} independent uniform \emph{p}-values has
#' distribution function \eqn{F_{k_i}(p) = 1 - (1 - p)^{k_i}}, so each
#' study \emph{p}-value is transformed by
#' \deqn{p_i \mapsto 1 - (1 - p_i)^{k_i}}
#' before combination. The adjusted
#' \emph{p}-values are again uniformly distributed on \eqn{[0, 1]} under the
#' null, so the null distribution used by the combination is unchanged.
#' 
#'
#' @return A numeric vector of combined \emph{p}-values corresponding 
#'   to each value of \code{mu}.
#'   
#' @importFrom stats pnorm
#'
#' @export
#'
#' @references
#' Tippett LHC. *Methods of Statistics*. Williams Norgate; 1931.
#' \doi{10.2307/3606890}
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
#' # Using Tippett's method to calculate the combined p-value
#' p_tippett(
#'      estimates = estimates,
#'      SEs = SEs,
#'      mu = mu,
#'      heterogeneity = "none",
#'      output_p = "two.sided",
#'      input_p = "greater"
#' )
#' 
#' Best-of-k adjustment: each study is the most significant of k = 3 experiments
#' p_tippett(
#'      estimates = estimates,
#'      SEs = SEs,
#'      mu = mu,
#'      input_p = "greater",
#'      k = rep(3, n)
#' )
p_tippett <- function(
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
  
  
  # Calculate the Tippett statistic
  ## sp <- 1 - (1 - apply(p, 2L, min))^n
  min_p <- apply(p, 2L, min)
  sp <- -expm1(n * log1p(-min_p))
  
  # Symmetrize to two-sided ONLY if requested and if inputs weren't already two-sided
  if (output_p == "two.sided" && input_p != "two.sided") {
    sp <- 2 * pmin(sp, 1 - sp)
  }
  return(sp)
}

