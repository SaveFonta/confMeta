#' Calculate the p-value
#'
#' @template thetahat
#' @template se 
#' @template w 
#' @param mu A numeric vector containing null hypothesis value(s).
#' @template phi
#' @template tau2
#' @template alternative 
#' @template distr
#' @template heterogeneity 
#' @param bound If \code{FALSE} (default), p-values that cannot be computed are reported as \code{NaN}.
#' If \code{TRUE}, they are reported as "> bound". This might, however, lead to situations in which
#' some of the p-values could be computed whereas others are reported as bounds. As a consequence, the
#' returned vector in such cases has class \code{"character"}. If the output vector must be of class \code{numeric}
#' use \code{bound = FALSE}.
#' @return Returns the p-value from the harmonic mean chi-squared test
#' based on study-specific estimates and standard errors.
#' @importFrom stats pchisq pf
#' @export
hMeanChiSqMu <- function(thetahat, se, 
                         w = rep(1, length(thetahat)),
                         mu = 0,
                         phi = estimatePhi(thetahat, se),
                         tau2 = estimateTau2(thetahat, se),
                         alternative = "none",
                         distr = c("chisq", "f"),
                         heterogeneity = c("additive", "multiplicative"),
                         bound = FALSE){
  
  # Check inputs
  stopifnot(is.numeric(thetahat),
            length(thetahat) > 0L,
            is.finite(thetahat),
            
            is.numeric(se),
            length(se) == 1L || length(se) == length(thetahat),
            is.finite(se),
            min(se) > 0,
            
            is.numeric(w),
            length(w) == length(thetahat),
            is.finite(w),
            min(w) > 0,
            
            is.numeric(phi) || is.null(phi),
            length(phi) == 1L,
            is.finite(phi) || is.null(phi),
            
            is.numeric(tau2) || is.null(tau2),
            length(tau2) == 1L,
            is.finite(tau2) || is.null(tau2),
            0 <= tau2 || is.null(tau2),
            
            !is.null(phi) || !is.null(tau2),
            
            is.numeric(mu),
            length(mu) > 0L,
            is.finite(mu),
            
            is.logical(bound),
            length(bound) == 1L,
            is.finite(bound),
            
            !is.null(alternative),
            length(alternative) == 1L,
            alternative %in% c("greater", "less", "two.sided", "none"),
            
            !is.null(distr),
            length(distr) == 1L,
            
            !is.null(heterogeneity),
            length(heterogeneity) == 1L)
  
  # match arguments
  distr <- match.arg(distr, several.ok = FALSE)
  if(length(se) == 1L) se <- rep(se, length(thetahat))
  
  # estimate heterogeneity
  if(heterogeneity == "additive"){
    denominator <- sqrt(se^2 + tau2)
  } else if(heterogeneity == "multiplicative"){
    denominator <- sqrt(se^2 * phi)
  }
  
  # store lengths of input vector
  n <- length(thetahat)
  m <- length(mu)
  
  # Case: Alternative is one of c("greater", "less", "two.sided")
  if(alternative != "none"){
    sw <- sum(sqrt(w))^2
    z <- vapply(seq_along(mu), function(i) (thetahat - mu[i]) / denominator, double(length(thetahat)))
    if(!is.matrix(z)) z <- as.matrix(z)
    zH2 <- vapply(seq_along(mu), function(i) sw / sum(w / z[, i]^2), double(1L))
    if(distr == "chisq"){
      res <- stats::pchisq(zH2, df = 1L, lower.tail = FALSE)
    } else if(distr == "f") {
      res <- stats::pf(zH2, df1 = 1L, df2 = (length(thetahat) - 1L), lower.tail = FALSE)
    }
    check_greater <- vapply(seq_len(ncol(z)), function(i) min(z[, i]) >= 0, logical(1L))
    check_less <- vapply(seq_len(ncol(z)), function(i) max(z[, i]) <= 0, logical(1L))
    break_p <- 1 / (2^n)
    if(alternative == "greater"){
      if(bound){
        res <- ifelse(check_greater, res / (2^n), paste(">", format(break_p, scientific = FALSE)))
        #res <- if(check_greater) res / (2^n) else paste(">", format(break_p, scientific = FALSE))
      } else {
        res <- ifelse(check_greater | check_less, res / (2^n), NaN)
        #res <- if(check_greater || check_less) res/(2^n) else NaN
      }
    }
    if(alternative == "less"){
      if(bound){
        res <- ifelse(check_less, res / (2^n), paste(">", format(break_p, scientific = FALSE)))
        #res <- if(check_less) res / (2^n) else paste(">", format(break_p, scientific = FALSE))
      } else {
        res <- ifelse(check_greater | check_less, res / (2^n),  NaN)
        #res <- if(check_greater || check_less) res / (2^n) else NaN
      }
    }
    if(alternative == "two.sided"){
      if(bound){
        res <- ifelse(check_greater | check_less, res / (2^(n - 1)), paste(">", format(2*break_p, scientific = FALSE)))
        #res <- if(check_greater || check_less) res / (2^(n - 1)) else paste(">", format(2*break_p, scientific = FALSE))
      } else {
        res <- ifelse(check_greater | check_less, res / (2^(n - 1)), NaN)
        #res <- if(check_greater || check_less) res / (2^(n - 1)) else NaN
      }
    }
  # Case when alternative is "none"
  } else if(alternative == "none") {
    sw <- sum(sqrt(w))^2
    zH2 <- vapply(seq_along(mu), function(i){
      z <- (thetahat - mu[i]) / denominator
      sw / sum(w / z^2)
    }, double(1L))
    if(distr == "chisq"){
      res <- stats::pchisq(q = zH2, df = 1, lower.tail = FALSE)
    } else if(distr == "f") {
      res <- stats::pf(zH2, df1 = 1, df2 = (length(thetahat) - 1), lower.tail = FALSE)
    }
  }
  return(res)
}
