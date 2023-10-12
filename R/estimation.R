#' @title Estimation of between-study heterogeneity
#' @rdname estimate_heterogeneity
#' @order 2
#'
#' @description \code{estimate_phi} estimates the the between-study heterogeneity
#'     \eqn{\phi}{phi} using a multiplicative model. The function is a
#'     modified version of Page 3 of: \cr
#'     Accounting for heterogeneity in meta-analysis using a
#'     multiplicative model - an empirical study,
#'     Mawdsley D. etal. 2016. DOI: 10.1002/jrsm.1216 \cr \cr
#'     The difference to the version in the paper is that the function,
#'     does not truncate the MSE at 1.
#'
# @note In the paper the weights should be squared!?
#'
#' @template estimates
#' @template SEs
#' @return For \code{estimate_phi}: the estimated heterogeneity parameter
#'     \eqn{\phi}{phi}, a \code{numeric} vector of length 1.
#'
#' @export
#'
#' @examples
#'     # Estimating heterogeneity using the multiplicative model
#'     estimate_phi(estimates = estimates, SEs = SEs)
estimate_phi <- function(estimates, SEs) {

    m <- stats::lm(estimates ~ 1, weights = 1 / SEs^2)
    mse <- stats::anova(m)$`Mean Sq`[1]
    ## corresponds to:
    ## mse <- sum(summary(m)$residuals^2) / (length(thetahat) - 1)
    ## max(1, mse) Ignore truncation at this point
    return(mse)
}


# Additive heterogeneity -------------------------------------------------------

#' @rdname estimate_heterogeneity
#' @order 1
#'
#' @description \code{estimate_tau2} estimates the between-study heterogeneity
#'     \eqn{\tau^{2}}{tau^2} using an additive model. The resulting parameter
#'     \eqn{\tau^2}{tau^2} is estimated through a call to
#'     \code{\link[meta]{metagen}} with \code{TE = estimates} and
#'     \code{seTE = SEs}. Other arguments to \code{\link[meta]{metagen}} can be
#'     passed via the \code{...} argument. If no arguments are passed via
#'     \code{...}, the following defaults are applied.
#'     \itemize{
#'       \item{\code{sm} = \code{"md"}}
#'       \item{\code{method.tau} = \code{"REML"}}
#'       \item{\code{control} = \code{list(maxiter = 1e5, stepadj = 0.25)}}
#'     }
#'     However, each of these defaults can be overwritten via \code{...}.
#'
#' @param ... Further arguments that are passed to \code{\link[meta]{metagen}}.
#'
#' @return For \code{estimate_tau2}: the estimated heterogeneity parameter
#'     \eqn{\tau^2}{tau^2}, a \code{numeric} vector of length 1.
#'
#' @export
#'
#' @examples
#'     # Example estimates and std. errors
#'     estimates <- c(0.21, 0.53, 0.24, 0.32, 0.19)
#'     SEs <- c(0.19, 0.39, 0.7, 1, 0.97)
#'
#'     # Estimating heterogeneity using the multiplicative model
#'     estimate_tau2(estimates = estimates, SEs = SEs)
estimate_tau2 <- function(estimates, SEs, ...) {

    # get names of ... arguments
    dotargs <- as.list(match.call()[-1])
    dotargs <- dotargs[!(names(dotargs) %in% c("estimates", "SEs"))]
    argnames <- names(dotargs)

    # set default arguments that can be overwritten by user via the ... argument
    # default args to be set
    default_args <- c(
        "sm" = "sm",
        "method.tau" = "method.tau",
        "control" = "control"
    )
    # keep only those that are not set by user
    default_args <- default_args[!(default_args %in% argnames)]
    # put all of the not overwritten defaults into a list
    args <- lapply(default_args, function(x) {
        if (x == "sm") value <- "MD"
        if (x == "method.tau") value <- "REML"
        if (x == "control") value <- list(maxiter = 1e5, stepadj = 0.25)
        return(value)
    })

    # set thetahat and se to the list
    args$TE <- estimates
    args$seTE <- SEs

    # append user set arguments
    args <- append(args, dotargs)

    # generate call to metagen and evaluate
    cc <- eval(as.call(append(list(meta::metagen), args)))

    return(cc$tau2)
}
