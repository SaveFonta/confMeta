#' Calculate the p-value using the harmonic mean chi-squared test.
#' 
#' @details
#' The function is vectorized over the argument \code{mu}.
#'
#' @template thetahat
#' @template se 
#' @template mu
#' @template phi
#' @template tau2
#' @template heterogeneity
#' @template alternative 
#' @template check_inputs
#' @template w 
#' @template distr
#' @return Returns the p-value from the harmonic mean chi-squared test
#' based on study-specific estimates and standard errors.
#' @export
hMeanChiSqMu <- function(
    thetahat,
    se,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    alternative = "none",
    check_inputs = TRUE,
    w = rep(1, length(thetahat)),
    distr = "chisq"
) {

  # Check inputs
  if (check_inputs) {
    check_inputs_p_value(
      thetahat = thetahat,
      se = se,
      mu = mu,
      heterogeneity = heterogeneity,
      alternative = alternative,
      phi = phi,
      tau2 = tau2
    )
    check_distr_arg(distr = distr)
    check_w_arg(w = w, thetahat = thetahat)
  }
  
  # match arguments
  if (length(se) == 1L) se <- rep(se, length(thetahat))
  
  # adjust se based on heterogeneity model
  se <- adjust_se(
    se = se,
    heterogeneity = heterogeneity,
    phi = phi,
    tau2 = tau2
  )
  
  # store lengths of input vector
  n <- length(thetahat)
  m <- length(mu)
  
  # Case: Alternative is one of c("greater", "less", "two.sided")
  if (alternative != "none") {
    sw <- sum(sqrt(w))^2
    z <- vapply(
      mu,
      function(mu) (thetahat - mu) / se,
      double(length(thetahat))
    )
    if (!is.matrix(z)) {
      z <- as.matrix(z)
    }
    zH2 <- vapply(
      seq_along(mu),
      function(i) sw / sum(w / z[, i]^2),
      double(1L)
    )
    if (distr == "chisq"){
      res <- stats::pchisq(zH2, df = 1L, lower.tail = FALSE)
    } else if (distr == "f") {
      res <- stats::pf(
        zH2,
        df1 = 1L,
        df2 = (length(thetahat) - 1L),
        lower.tail = FALSE
      )
    }
    check_greater <- vapply(
      seq_len(ncol(z)),
      function(i) min(z[, i]) >= 0,
      logical(1L)
    )
    check_less <- vapply(
      seq_len(ncol(z)),
      function(i) max(z[, i]) <= 0,
      logical(1L)
    )
    break_p <- 1 / (2^n)
    if (alternative == "greater") {
      res <- ifelse(check_greater | check_less, res / (2^n), NaN)
    }
    if (alternative == "less") {
      res <- ifelse(check_greater | check_less, res / (2^n),  NaN)
    }
    if (alternative == "two.sided") {
      res <- ifelse(check_greater | check_less, res / (2^(n - 1)), NaN)
    }
  # Case when alternative is "none"
  } else if (alternative == "none") {
    sw <- sum(sqrt(w))^2
    zH2 <- vapply(
      mu,
      function(mu){
        z <- (thetahat - mu) / se
        sw / sum(w / z^2)
      },
      double(1L)
    )
    res <- if (distr == "chisq") {
      stats::pchisq(q = zH2, df = 1, lower.tail = FALSE)
    } else if (distr == "f") {
      stats::pf(zH2, df1 = 1, df2 = (length(thetahat) - 1), lower.tail = FALSE)
    }
  }
  return(res)
}