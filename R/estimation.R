# Multiplicative heterogeneity -------------------------------------------------------------------------------


#' Estimate phi for multiplicative heterogeneity
#'
#' @description 
#' A modified version of Page 3 of: \cr \cr
#' Accounting for heterogeneity in meta-analysis using a multiplicative model -
#' an empirical study, Mawdsley D. etal. 2016. DOI: 10.1002/jrsm.1216
#' 
#' This function, in contrast to the paper, does not truncate the MSE at 1.
#'
#' @note In the paper the weights should be squared!?
#' 
#' @template thetahat 
#' @template se
#' @return Estimated phi, a \code{numeric} vector of length 1.
#' 
#' @importFrom stats lm anova
#' 
#' @export
#' 
#' @examples
#' thetahat <- c(0.21, 0.53, 0.24, 0.32, 0.19)
#' se <- c(0.19, 0.39, 0.7, 1, 0.97)
#' estimatePhi(thetahat = thetahat, se = se)
#'
#' thetahat[1] <- 7
#' estimatePhi(thetahat = thetahat, se = se)
estimatePhi <- function(thetahat, se){
  
  m <- stats::lm(thetahat ~ 1, weights = 1 / se^2)
  mse <- stats::anova(m)$`Mean Sq`[1]
  ## corresponds to:
  ## mse <- sum(summary(m)$residuals^2) / (length(thetahat) - 1)
  ## max(1, mse) Ignore truncation at this point
  return(mse)
}


# Additive heterogeneity -------------------------------------------------------------------------------

#' Estimate tau2 for additive heterogeneity
#' 
#' @description
#' This function estimates tau2 for the additive heterogeneity model. The estimation is done via a call to the function \code{\link[meta]{metagen}} which offers a lot of
#' flexibility in terms of how tau2 should be estimated. Arguments can be passed via the \code{...} argument. Note, however, that the following arguments that are passed 
#' to the \code{\link[meta]{metagen}} function are set within this function: 
#' \itemize{
#'   \item{\code{sm} = \code{"md"}}
#'   \item{\code{method.tau} = \code{"REML"}}
#'   \item{\code{control} = \code{list(maxiter = 1e5, stepadj = 0.25)}}
#' }
#' 
#' However, all of these arguments can be overwritten by the user via the \code{...} argument.
#' 
#' @template thetahat
#' @template  se
#' @param ... Further arguments that are passed to \code{\link[meta]{metagen}}.
#' @return Estimated tau2, a \code{numeric} vector of length 1.
#' 
#' @importFrom meta metagen
#' 
#' @export
#' 
#' @examples
#' thetahat <- c(0.65, -0.5, 0.03, -0.51, 0.85, 0.75, -0.22, 1.02, 0.36, 0.44)
#' se <- c(0.2, 0.2, 0.19, 0.22, 0.21, 0.18, 0.17, 0.18, 0.21, 0.19)
#' estimateTau2(thetahat = thetahat, se = se)
#'
#' thetahat[1] <- 7
#' estimateTau2(thetahat = thetahat, se = se)
estimateTau2 <- function(thetahat, se, ...){
  
  # get names of ... arguments
  dotargs <- as.list(match.call()[-1])
  dotargs <- dotargs[!(names(dotargs) %in% c("thetahat", "se"))]
  argnames <- names(dotargs)
  
  # set default arguments that can be overwritten by user via the ... argument
  # default args to be set
  default_args <- c("sm" = "sm", 
                    "method.tau" = "method.tau",
                    "control" = "control")
  default_args <- default_args[!(default_args %in% argnames)] # keep only those that are not set by user
  args <- lapply(default_args, function(x){ # put all of the not overwritten defaults into a list
    if(x == "sm") value <- "MD"
    if(x == "method.tau") value <- "REML"
    if(x == "control") value <- list(maxiter = 1e5, stepadj = 0.25)
    return(value)
  })
  
  # set thetahat and se to the list
  args$TE <- thetahat
  args$seTE <- se
  
  # append user set arguments
  args <- append(args, dotargs)
  
  # generate call to metagen and evaluate
  cc <- eval(as.call(append(list(meta::metagen), args)))
  
  return(cc$tau2)
}
