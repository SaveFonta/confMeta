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
#' @template check_inputs
#' @return Returns the p-value from the harmonic mean chi-squared test
#' based on study-specific estimates and standard errors.
#' @importFrom stats pchisq pf
#' @export
hMeanChiSqMu <- function(thetahat, se, 
                         w = rep(1, length(thetahat)),
                         mu = 0,
                         phi = NULL,
                         tau2 = NULL,
                         alternative = "none",
                         distr = "chisq",
                         heterogeneity = "none",
                         check_inputs = TRUE){
  # Check inputs
  if(check_inputs){
    stopifnot(is.numeric(thetahat),
              length(thetahat) > 0L,
              is.finite(thetahat),
              
              is.numeric(se),
              is.finite(se),
              min(se) > 0,
              
              is.numeric(w),
              length(w) == length(thetahat),
              is.finite(w),
              min(w) > 0,
              
              is.numeric(mu),
              length(mu) > 0L,
              is.finite(mu),
              
              !is.null(alternative),
              length(alternative) == 1L,
              alternative %in% c("greater", "less", "two.sided", "none"),
              
              !is.null(distr),
              length(distr) == 1L,
              
              !is.null(heterogeneity),
              length(heterogeneity) == 1L)
    
    if(!length(se) %in% c(1L, length(thetahat))) stop("Length of argument 'se' must be either 1 or length(thetahat).")
    if(!is.null(phi) && length(phi) > 1L) stop("Argument 'phi' must be of length 1.")
    if(!is.null(tau2) && length(tau2) > 1L) stop("Argument 'tau2' must be of length 1.")
    if(is.null(phi) && heterogeneity == "multiplicative") stop("If heterogeneity = 'multiplicative', phi must be provided.")
    if(is.null(tau2) && heterogeneity == "additive") stop("If heterogeneity = 'additive', tau2 must be provided.")
    if(!distr %in% c("f", "chisq")) stop("Argument 'distr' must be one of c('f', 'chisq').")
    if(!is.null(phi) && length(phi) != 1L) stop("Argument 'phi' must be of length 1.")
    if(!is.null(tau2) && length(tau2) != 1L) stop("Argument 'tau2' must be of length 1.")
    if(!is.null(phi) && heterogeneity == "additive") warning("If heterogeneity = 'additive', argument 'phi' is ignored.")
    if(!is.null(tau2) && heterogeneity == "multiplicative") warning("If heterogeneity = 'multiplicative', argument 'tau2' is ignored.")
    if(heterogeneity == "none" && (!is.null(phi) || !is.null(tau2))) warning("If heterogeneity = 'none', arguments 'tau2' and 'phi' are ignored.")
    if(!is.null(phi) && (!is.finite(phi) || phi < 0)) stop("Argument 'phi' must be finite and larger than 0.")
  }
  
  # match arguments
  if(length(se) == 1L) se <- rep(se, length(thetahat))
  
  # add heterogeneity
  if(heterogeneity == "additive") denominator <- sqrt(se^2 + tau2)
  if(heterogeneity == "multiplicative") denominator <- sqrt(se^2 * phi)
  if(heterogeneity == "none") denominator <- se
  
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
      res <- ifelse(check_greater | check_less, res / (2^n), NaN)
    }
    if(alternative == "less"){
      res <- ifelse(check_greater | check_less, res / (2^n),  NaN)
    }
    if(alternative == "two.sided"){
      res <- ifelse(check_greater | check_less, res / (2^(n - 1)), NaN)
    }
  # Case when alternative is "none"
  } else if(alternative == "none") {
    sw <- sum(sqrt(w))^2
    zH2 <- vapply(seq_along(mu), function(i){
      z <- (thetahat - mu[i]) / denominator
      sw / sum(w / z^2)
    }, double(1L))
    res <- if(distr == "chisq"){
      stats::pchisq(q = zH2, df = 1, lower.tail = FALSE)
    } else if(distr == "f") {
      stats::pf(zH2, df1 = 1, df2 = (length(thetahat) - 1), lower.tail = FALSE)
    }
  }
  return(res)
}
