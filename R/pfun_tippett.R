#' @title Tippett's method
#' @family p-value combination functions
#'
#' @description
#' Tippett's method for combining \emph{p}-values across studies. 
#' This method evaluates the minimum individual study \emph{p}-value.
#'
#' @param estimates Numeric vector of study-level effect estimates.
#' @param SEs Numeric vector of corresponding standard errors.
#' @param mu Numeric scalar or vector of null values for the overall effect
#'     (default: 0). 
#' @param heterogeneity Character string: \code{"none"} (default),
#'     \code{"additive"}, or \code{"multiplicative"}. Determines whether
#'     standard errors are adjusted for between-study heterogeneity using
#'     \code{tau2} or \code{phi}.
#' @param phi Multiplicative heterogeneity parameter (if applicable).
#' @param tau2 Additive heterogeneity parameter (if applicable).
#' @param check_inputs Logical (default \code{TRUE}). If \code{TRUE},
#'     perform input validation.
#' @param output_p Character string specifying the combined
#'     \emph{p}-value type: \code{"two.sided"} (default) or \code{"one.sided"}.
#'     This controls whether the final combined \emph{p}-value is symmetrized.
#'     \strong{Note:} To construct valid \emph{p}-value functions, the \code{\link{confMeta}}
#'     function strictly requires \code{"two.sided"}.
#' @param input_p Type of study-level \emph{p}-values used in the combination:
#'     \code{"greater"} (default), \code{"less"}, or \code{"two.sided"}.
#'     If \code{"greater"} or \code{"less"}, one-sided \emph{p}-values are
#'     combined.
#'
#' @details
#' The Tippett combined \emph{p}-value, \eqn{p_T}, for \eqn{k} studies is 
#' defined directly as:
#' \deqn{p_T = 1 - (1 - \min\{p_1, \dots, p_k\})^k}
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
p_tippett <- function(
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
  
  # Calculate the Tippett statistic
  sp <- 1 - (1 - apply(p, 2L, min))^n
  
  # Symmetrize to two-sided ONLY if requested and if inputs weren't already two-sided
  if (output_p == "two.sided" && input_p != "two.sided") {
    sp <- 2 * pmin(sp, 1 - sp)
  }
  return(sp)
}
