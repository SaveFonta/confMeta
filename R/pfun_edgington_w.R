#' @title Weighted Edgington’s method
#' @rdname p_value_functions
#' @order 8
#'
#' @description
#' Weighted generalization of Edgington’s method for combining
#' \emph{p}-values across studies. The method forms a weighted sum
#' of individual study \emph{p}-values and evaluates it against the
#' exact or approximate null distribution.
#'
#' If all weights are equal, this reduces to the classical Edgington
#' procedure, where the null distribution is given by the Irwin–Hall law.
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
#' study \emph{p}-values.  
#'
#' Under the null hypothesis, the distribution of \eqn{S} can be obtained
#' in two ways:
#' \itemize{
#'   \item For small numbers of studies or unbalanced weights, the exact
#'   distribution is computed via the Barrow–Smith inclusion–exclusion formula.
#'   \item For sufficiently many studies (or balanced weights), the distribution
#'   is approximated by a normal distribution with mean \eqn{\sum w_i / 2}
#'   and variance \eqn{\sum w_i^2 / 12}.
#' }
#'
#' The condition for validity of the normal approximation can be expressed
#' in terms of the **effective sample size**:
#' \deqn{n_{\mathrm{eff}} = \frac{\|w\|_2^4}{\|w\|_4^4}.}
#' If \eqn{n_{\mathrm{eff}} \geq 12}, the approximation error is small
#'  When all weights are equal, this reduces to
#' the simple rule \eqn{n \geq 12}.
#'
#' By construction, the output is always two-sided:
#' \deqn{p_{2s} = 2 \min(p, 1-p).}
#'
#' @return A numeric vector of combined \emph{p}-values corresponding
#'     to each value of \code{mu}.


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
    confMeta:::check_inputs_p_value(estimates = estimates, SEs = SEs, 
                                    mu = mu, heterogeneity = heterogeneity, 
                                    phi = phi, tau2 = tau2)
    confMeta:::check_alternative_arg_edg(alternative = alternative)
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
    SEs <- confMeta:::adjust_se(SEs = SEs, heterogeneity = heterogeneity,
                                phi = phi, tau2 = tau2)
  }
  
  # Compute z-values (matrix: n x length(mu))
  z <- confMeta:::get_z(estimates = estimates, SEs = SEs, mu = mu)
  
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
    sp <- confMeta:::pirwinhall(q = colSums(p), n = n, approx = approx)
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
      if (!use_normal && n > 18) {
        stop("Exact method infeasible for n > 18; use approx=TRUE or lower weights threshold") #cause with 2^18 R can crush
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

# Vectorize over mu
#p_edgington_w2 <- Vectorize(FUN = p_edgington_w2., vectorize.args = "mu")

"completelY useless, why did he do it?"



# NEXT STEP: write a testhat

