#' @param heterogeneity One of \code{c("none", "additive", "multiplicative")}. If \code{heterogeneity = "none"} p-values are returned for the 
#' passed \code{se} without any adaption. If \code{heterogeneity = "additive"}, the standard errors \code{se} are reassigned the value of \code{sqrt(se^2 + tau2)} before
#' computation of the p-values. If \code{heterogeneity = "multiplicative"}, the standard errors \code{se} are multiplied with the value of \code{phi} before
#' computation of the p-values. Defaults to \code{"none"}.
