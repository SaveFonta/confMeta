#' Implementation of the "k-Trials" rule.
#' @details
#' The function is vectorized over the argument \code{mu}.
#' 
#' @template thetahat
#' @template se
#' @template mu
#' @template tau2
#' @template phi
#' @template heterogeneity
#' @template check_inputs
#'
#' @return The corresponding p-value given mu under the null-hypothesis.
#' @export
#'
#' @examples
#' thetahat <- c(0.38041824, -0.22170681, -0.09155859)
#' se <- c(0.2256202, 0.2432796, 0.1051743)
#' phi <- estimatePhi(thetahat = thetahat, se = se)
#' tau2 <- estimateTau2(
#'     thetahat,
#'     se,
#'     control = list(stepadj = 0.5, maxiter = 1000, threshold = 1e-6)
#' )
#' mymu <- seq(
#'     min(thetahat - 3 * se),
#'     max(thetahat + 3 * se),
#'     length.out = 1000
#' )
#' p_add <- kTRMu(
#'     thetahat, se,
#'     mymu,
#'     tau2 = tau2,
#'     phi = phi,
#'     heterogeneity = "additive"
#' )
#' p_mult <- kTRMu(
#'     thetahat,
#'     se,
#'     mymu,
#'     tau2 = tau2,
#'     phi = phi,
#'     heterogeneity = "multiplicative"
#' )
kTRMu <- function(
  thetahat,
  se,
  mu,
  phi = NULL,
  tau2 = NULL,
  heterogeneity = "none",
  check_inputs = TRUE
) {

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
  }
  
  # recycle `se` if needed
  if (length(se) == 1L) se <- rep(se, length(thetahat))
  
  # adjust se based on heterogeneity model
  se <- adjust_se(
    se = se,
    heterogeneity = heterogeneity,
    phi = phi,
    tau2 = tau2
  )
  
  # Get lengths
  n <- length(thetahat)
  
  if(alternative == "none") {
    z <- vapply(mu, function(mu) (thetahat - mu) / se, double(length(thetahat)))
    if (is.null(dim(z)))
      dim(z) <- c(1L, length(z))
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE) # ReplicationSuccess::z2p
    pmax <- apply(p, 2L, max)
  } else {
    stop("Invalid argument 'alternative'.")
  }
  
  return(pmax^n)
}
