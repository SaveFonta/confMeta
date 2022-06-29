#' Implementation of the "k-Trials" rule.
#'
#' @template thetahat 
#' @template se 
#' @param mu A numeric vector containing null hypothesis value(s).
#' @template tau2 
#' @template phi 
#' @template heterogeneity
#'
#' @return The corresponding p-value given mu under the null-hypothesis
#' @importFrom ReplicationSuccess z2p
#' @export
#'
#' @examples
#' thetahat <- c(0.38041824, -0.22170681, -0.09155859)
#' se <- c(0.2256202, 0.2432796, 0.1051743)
#' phi <- estimatePhi(thetahat = thetahat, se=se)
#' tau2 <- estimateTau2(thetahat, se, control = list(stepadj = 0.5, maxiter = 1000, threshold = 1e-6))
#' mymu <- seq(min(thetahat - 3 * se), max(thetahat + 3 * se), length.out=1000)
#' p_add <- kTRMu(thetahat, se, mymu, tau2 = tau2, phi = phi, heterogeneity = "additive")
#' p_mult <- kTRMu(thetahat, se, mymu, tau2 = tau2, phi = phi, heterogeneity = "multiplicative")
kTRMu <- function(thetahat, se, 
                  mu, 
                  phi = NULL, 
                  tau2 = NULL,
                  heterogeneity = c("none", "additive", "multiplicative")){
  
  if(is.null(tau2) && heterogeneity == "additive") stop("If heterogeneity = 'additive', tau2 must be provided")
  if(is.null(phi) && heterogeneity == "multiplicative") stop("If heterogeneity = 'multiplicative', phi must be provided")
  if(heterogeneity == "none" && (!is.null(tau2) || !is.null(phi))) warning(paste0("Ignoring parameter(s) phi and/or tau2 as heterogeneity = 'none'"))
  
  if(heterogeneity == "additive")
    se <- sqrt(se^2+tau2)
  if(heterogeneity == "multiplicative")
    se <- se*sqrt(phi)
  z <- vapply(mu, function(mu) (thetahat - mu)/se, double(length(thetahat)))
  if(is.null(dim(z)))
    dim(z) <- c(1, length(z))
  p <- ReplicationSuccess::z2p(z, alternative="two.sided")
  dim(p) <- dim(z)
  pmax <- apply(p, 2, max)
  k <- length(thetahat)
  return(pmax^k)
}
