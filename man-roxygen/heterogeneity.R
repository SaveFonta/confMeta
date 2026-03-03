#' @param heterogeneity One of \code{c("none", "additive", "multiplicative")}. 
#'     If \code{heterogeneity = "none"}, \emph{p}-values are returned for the 
#'     passed \code{SEs} without any adaptation. 
#'     If \code{heterogeneity = "additive"}, the standard errors are reassigned 
#'     the value \eqn{\sqrt{SEs^2 + \text{tau2}}} before computation of the 
#'     \emph{p}-values. 
#'     If \code{heterogeneity = "multiplicative"}, the standard errors \code{SEs} 
#'     are multiplied by \eqn{\sqrt{\text{phi}}} before computation of the 
#'     \emph{p}-values. 
#'     Defaults to \code{"none"}.
