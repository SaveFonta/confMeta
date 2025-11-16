#' @title Weighted Edgington’s method
#' @rdname p_value_functions
#' @order 8
#'
#'@export
#'
#' @description
#' Weighted generalization of Edgington’s method for combining
#' \emph{p}-values across studies. The method forms a weighted sum
#' of individual study \emph{p}-values and evaluates it against the
#' exact or approximate null distribution.
#'
#' If all weights are equal, this reduces to the classical Edgington
#' procedure, where the null distribution is given by the Irwin–Hall distribution.
#'
#' @param estimates Numeric vector of study-level effect estimates.
#' @param SEs Numeric vector of corresponding standard errors.
#' @param mu Numeric scalar or vector of null values for the overall effect
#'     (default: 0). The function is vectorized over \code{mu}.
#' @param alternative Either \code{"two.sided"} (default), or \code{"one.sided"}.
#'     **Currently has no effect on the output**: all results are converted
#'     to two-sided combined \emph{p}-values, regardless of the argument.
#'     The argument is included only for consistency with other functions (since this function is a blueprint of p_edgington)
#'     and for input validation.
#' @param approx Logical (default \code{TRUE}). If \code{TRUE}, use a normal
#'     approximation for the weighted sum when the condition defined by \code{approx_rule}
#'     is met
#' @param input_p Type of study-level \emph{p}-values used in the combination:
#'     \code{"greater"}, \code{"less"}, or \code{"two.sided"}.
#'     If \code{"greater"} or \code{"less"}, one-sided \emph{p}-values are
#'     first computed, then symmetrized to two-sided in the output.
#' @param w Numeric vector of nonnegative weights, same length as \code{estimates}.
#'     Defaults to equal weights.
#' @param heterogeneity Character string: \code{"none"} (default),
#'     \code{"additive"}, or \code{"multiplicative"}. Determines whether
#'     standard errors are adjusted for between-study heterogeneity using
#'     \code{tau2} or \code{phi}.
#' @param phi Multiplicative heterogeneity parameter (if applicable).
#' @param tau2 Additive heterogeneity parameter (if applicable).
#' @param approx_rule Rule for normal approximation: \code{"n"} 
#'     uses the number of studies; \code{"neff"} (default) uses the effective sample
#'     size criterion (see Details).
#' @param neff_cut Numeric threshold (default 12). If \code{approx_rule="n"},
#'     normal approximation is used when \eqn{n \geq 12}. If
#'     \code{approx_rule="neff"}, normal approximation is used when
#'     \eqn{\|w\|_2^4 / \|w\|_4^4 \geq 12}.
#' @param check_inputs Logical (default \code{TRUE}). If \code{TRUE},
#'     perform input validation.
#'
#' @details
#' The weighted Edgington statistic is defined as
#' \deqn{S = \sum_{i=1}^k w_i p_i,}
#' where \eqn{w_i} are positive study weights and \eqn{p_i} are individual
#' study \emph{p}-values. Under the global null hypothesis, each \eqn{p_i}
#' is assumed to follow a \eqn{Unif(0, 1)} distribution. 
#' 
#' @section Null Distribution and Approximation:
#' The CDF of the test statistic S under the null, \eqn{F(t) = P(S \leq t)},
#' is computed in one of two ways:
#' \itemize{
#'   \item **Exact Method:** The function uses the exact Barrow-Smith
#'     inclusion-exclusion formula to compute the CDF.
#'     This method is computationally intensive and is infeasible for `n > 18`
#'     studies, at which point the function will stop with an error if
#'     `approx = FALSE` is used.
#'
#'   \item **Normal Approximation:** For a large number of studies or
#'     sufficiently balanced weights, S is approximated by a Normal
#'     distribution  with:
#'     \deqn{E[S] = \frac{1}{2}\sum_{i=1}^k w_i}
#'     \deqn{Var(S) = \frac{1}{12}\sum_{i=1}^k w_i^2}
#' }
#'
#'
#' @section Approximation Rule:
#' The `approx = TRUE` argument enables the normal approximation, but it is
#' only used if a condition (controlled by \code{approx_rule}) is met.
#' \itemize{
#'   \item \code{approx_rule = "n"}: Uses the approximation if the
#'     number of studies n > neff_cut.
#'   \item \code{approx_rule = "neff"} (default): Uses the approximation if the
#'     **effective sample size** n_eff > neff_cut.
#' }
#' The effective sample size is defined as:
#' \deqn{n_{\mathrm{eff}} = \frac{(\sum w_i^2)^2}{\sum w_i^4} = \frac{\|w\|_2^4}{\|w\|_4^4}}
#' The default threshold \code{neff_cut = 12} is based on error bounds, which indicate the approximation is
#' sufficiently accurate when this condition is met. This
#' rule is more robust than \code{approx_rule = "n"} when weights are
#' unbalanced.
#'
#'
#'
#' @section Approximation Error:
#' The \code{neff_cut} parameter directly controls the tolerance for
#' the approximation error.
#'   **Edgeworth Approximation:** This provides an
#'   *approximation* of the error, which directly relates to the
#'   \eqn{n_{\mathrm{eff}}} criterion used in this function:
#'   \deqn{\sup_{t \in \mathbb{R}} |F_n^w(t) - \Phi(t)| \approx \frac{||\phi^{(3)}||_{\infty}}{20}\frac{||w||_{4}^{4}}{||w||_{2}^{4}} \approx \frac{0.028}{n_{\mathrm{eff}}}}
#'
#' The default \code{neff_cut = 12} is chosen based on this, as it
#' corresponds to an approximate maximum error of
#' \eqn{0.028 / 12 \approx 0.0023}.
#' If you change \code{neff_cut}, you are changing this error tolerance.
#' 
#' 
#'
#' @section Output p-value:
#' The final output depends on the \code{input_p} argument:
#' \itemize{
#'   \item If \code{input_p} is \code{"greater"} or \code{"less"}, the
#'     input \eqn{p_i} are one-sided. The function computes the
#'     one-sided \emph{p}-value \eqn{sp = P(S \leq s_{\mathrm{obs}})} and
#'     returns a **symmetrized, two-sided \emph{p}-value**:
#'     \eqn{p_{2s} = 2 \min(sp, 1-sp)}.
#'   \item If \code{input_p} is \code{"two.sided"}, the input \eqn{p_i}
#'     are two-sided. The function returns \eqn{sp = P(S \leq s_{\mathrm{obs}})}
#'     **directly, without symmetrization**. Note that summing two-sided
#'     \emph{p}-values is not a standard application of Edgington's method.
#' }
#'
#' @return A numeric vector of combined \emph{p}-values corresponding
#'   to each value of \code{mu}. Note that the output is only
#'   symmetrized to a two-sided \emph{p}-value if \code{input_p} is
#'   \code{"greater"} or \code{"less"}.
#'   
#'   
#'   
#'   
#'   
#'    @references
#' D. L. Barrow and P. W. Smith. Spline notation applied to a volume
#' problem. *The American Mathematical Monthly*, 86(1):50-51, 1979.
#'
#' H. Cramér. *Mathematical Methods of Statistics*. Princeton University
#' Press, 1946.
#'
#' E. G. Olds. A note on the convolution of uniform distributions.
#' *The Annals of Mathematical Statistics*, 23(2):282-285, 1952.
#'
#' B. D. Ripley. *Stochastic Simulation*. John Wiley & Sons, Hoboken,
#' NJ, 1987.
#' 
#' Add others on p value functions







p_edgington_w <- function(
    estimates, SEs, mu = 0,
    alternative = "two.sided",        
    approx = TRUE,
    input_p = "greater",              
    w = rep(1, length(estimates)),
    heterogeneity = "none",           
    phi = NULL, tau2 = NULL,
    approx_rule = "neff",                
    neff_cut = 12,
    check_inputs = TRUE
){
  # -------------------------
  # Input checks
  # -------------------------
  if (check_inputs) {
    check_inputs_p_value(estimates = estimates, SEs = SEs, 
                                    mu = mu, heterogeneity = heterogeneity, 
                                    phi = phi, tau2 = tau2)
    check_alternative_arg_edg(alternative = alternative)
  }
  
  # Recycle SEs if only one provided
  if (length(SEs) == 1L) {
    SEs <- rep(SEs, length(estimates))
  }
  
  n <- length(estimates)
  
  stopifnot(length(SEs) == n)
  stopifnot(length(w) == n, all(is.finite(w)), all(w > 0), all(SEs > 0))
  
  # Optionally adjust SEs for heterogeneity
  if (heterogeneity != "none") {
    SEs <- adjust_se(SEs = SEs, heterogeneity = heterogeneity,
                                phi = phi, tau2 = tau2)
  }
  
  # Compute z-values (matrix: n x length(mu))
  z <- get_z(estimates = estimates, SEs = SEs, mu = mu)
  
  # Convert to p-values depending on input_p
  p <- switch(input_p,
              "two.sided" = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
              "greater"   = stats::pnorm(z, lower.tail = FALSE),
              "less"      = stats::pnorm(z, lower.tail = TRUE),
              stop("input_p must be 'greater','less','two.sided'")
  )
  
  #p is a matrix with --> each row a different study
  #each column p value associated to a different mu

    # Case 1: unweighted → Irwin–Hall distribution
  if (all(w == 1)) {
    sp <- pirwinhall(q = colSums(p), n = n, approx = approx)
  } else {
    # Case 2: weighted version
    # Decide whether to use normal approximation
    use_normal <- FALSE
    if (approx) {
      if (approx_rule == "n") {
        use_normal <- (n >= neff_cut)
      } else if (approx_rule == "neff") {
        s2 <- sum(w^2); s4 <- sum(w^4)
        neff <- (s2^2) / s4 
        use_normal <- (neff >= neff_cut)
      }
      else stop("approx_rule must be 'n' or 'neff'")
      
    }
    
    if (use_normal) {
      # Normal approximation for large n or neff
      mn <- sum(w) / 2
      sd <- sqrt(sum(w^2) / 12)
      sp <- stats::pnorm(q = as.numeric(colSums(p * w)), #weighted sum of p values for each column (i.e. each mu value), then take Normal cdf
                         mean = mn, sd = sd)
    } else {
      if (n > 18) {
        # We cannot do exact for n > 18
        
        if (!approx) {
          # 1) approx = FALSE: explicitly ask to use approx = TRUE
          stop(
            paste0(
              "Exact method infeasible for n > 18 in p_edgington_w with approx = FALSE.\n",
              "Please call the function again with `approx = TRUE`."
            )
          )
        } else {
          # 2) approx = TRUE but the approximation threshold was not met
          #    report weight ratio (and neff if available)
          w_ratio <- max(w) / min(w)
          
          msg <- paste0(
            "Exact method infeasible for n > 18 in weighted Edgington.\n",
            "The normal approximation was not used \n"
          )
          
          if (!is.na(neff)) {
            msg <- paste0(
              msg,
              "Effective sample size neff = ", signif(neff, 4),
              " < neff_cut = ", neff_cut, ".\n"
            )
          }
          
          msg <- paste0(
            msg,
            "Weight imbalance: max(w)/min(w) = ", signif(w_ratio, 4), ".\n",
            "Consider lowering the weight imbalance (e.g., modifying or rescaling the weights) ",
            "or relaxing the approximation threshold (e.g., a smaller `neff_cut`)."
          )
          
          stop(msg)
        }
      }
      # Exact calculation using Barrow–Smith inclusion–exclusion formula
      vertices <- as.matrix(expand.grid(rep(list(0:1), n))) #generate all vertices of the n dim hypercube (2^n vectors, where n =#studies)
      signv    <- (-1) ^ rowSums(vertices)
      wTv      <- as.numeric(vertices %*% w) #it is w'v where v are the vertices generated
      denom    <- factorial(n) * prod(w)
      
      cdf_weighted_uniform <- function(ew){
        if (ew <= 0) return(0)  #since P(T<=0) = 0 cause T is always positive since it is a sum of weighted U(0,1) w positive weights
        S <- sum(w) # since the weigthed sum T has support [0, sum(weights)], this mean sum of weights is the max possible value
        if (ew >= S) return(1) #since P(T<=S) = 1
        idx <- which(wTv <= ew) #identify which sum respect the condition of the Indicator in the formula
        if (length(idx) == 0) return(0)
        dif <- ew - wTv[idx] 
        sum(signv[idx] * dif^n) / denom
      }
      
      ew_all <- as.numeric(colSums(p * w))
      sp <- vapply(ew_all, cdf_weighted_uniform, numeric(1))
    }
  }
  
  # Final output: by default, convert to two-sided p-values
  # If input_p == "two.sided", p_i are already symmetric and no change is applied
  if (input_p != "two.sided") {
    sp <- 2*pmin(sp, 1 - sp)
  }
  return(sp)
}



