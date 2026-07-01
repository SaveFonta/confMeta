#' @title Estimation of between-study heterogeneity
#' @rdname estimate_heterogeneity
#' @order 2
#'
#' @description \code{estimate_phi} estimates the between-study heterogeneity
#'     \eqn{\phi} using a multiplicative model. The function is a
#'     modified version of Page 3 of: \cr
#'     Accounting for heterogeneity in meta-analysis using a
#'     multiplicative model -- an empirical study,
#'     Mawdsley D. et al. 2016. \doi{10.1002/jrsm.1216} \cr \cr
#'
#'
#' @template estimates
#' @template SEs
#' @return For \code{estimate_phi}: the estimated heterogeneity parameter
#'     \eqn{\phi}, a \code{numeric} vector of length 1.
#'
#' @export
#'
#' @examples
#'     # Example estimates and std. errors
#'     estimates <- c(0.21, 0.53, 0.24, 0.32, 0.19)
#'     SEs <- c(0.19, 0.39, 0.7, 1, 0.97)
#'     
#'     estimate_phi(estimates = estimates, SEs = SEs)
#' @importFrom stats lm anova
estimate_phi <- function(estimates, SEs) {

    m <- stats::lm(estimates ~ 1, weights = 1 / SEs^2)
    mse <- stats::anova(m)$`Mean Sq`[1]
    ## corresponds to:
    ## mse <- sum(summary(m)$residuals^2) / (length(thetahat) - 1)
    phi <- max(c(1, mse)) ## truncate at 1
    return(phi)
}


# Additive heterogeneity -------------------------------------------------------

#' @rdname estimate_heterogeneity
#' @order 1
#'
#' @description \code{estimate_tau2} estimates the between-study heterogeneity
#'     \eqn{\tau^{2}} using an additive model. The resulting parameter
#'     \eqn{\tau^2} is estimated through a call to
#'     \code{meta::metagen} with \code{TE = estimates} and
#'     \code{seTE = SEs}. Other arguments to \code{meta::metagen} can be
#'     passed via the \code{...} argument. If no arguments are passed via
#'     \code{...}, the following defaults are applied.
#'     \itemize{
#'       \item \code{sm = "MD"} 
#'       \item \code{method.tau = "REML"} 
#'       \item \code{control = list(maxiter = 1e5, stepadj = 0.25)} 
#'     }
#'     However, each of these defaults can be overwritten via \code{...}.
#'
#' @param ... Further arguments that are passed to \code{meta::metagen()}.
#'
#' @return For \code{estimate_tau2}: the estimated heterogeneity parameter
#'     \eqn{\tau^2}, a \code{numeric} vector of length 1.
#'
#' @export
#'
#' @examples
#'     # Example estimates and std. errors
#'     estimates <- c(0.21, 0.53, 0.24, 0.32, 0.19)
#'     SEs <- c(0.19, 0.39, 0.7, 1, 0.97)
#'
#'     # Estimating heterogeneity using the additive model
#'     estimate_tau2(estimates = estimates, SEs = SEs)
#' @importFrom meta metagen
estimate_tau2 <- function(estimates, SEs, ...) {
  
  # Capture dot arguments 
  args <- list(TE = estimates,
               seTE = SEs, ...)
  
  # Set defaults ONLY if the user hasn't provided them in ...
  if (is.null(args$sm)) args$sm <- "MD"
  if (is.null(args$method.tau)) args$method.tau <- "REML"
  if (is.null(args$control)) args$control <- list(maxiter = 1e5, stepadj = 0.25)
  
  
  # Run meta::metagen using do.call
  cc <- do.call(meta::metagen, args)
  
  return(cc$tau2)
}


# this is old code, maybe there is a smarter way to write it TO DO
