#' Calculate confidence intervals based on the harmonic mean chi-squared test
#' 
#' @template thetahat
#' @template se
#' @template w
#' @template phi
#' @template tau2
#' @param level Numeric vector of length 1 specifying the level of the confidence interval. Defaults to 0.95.
#' @template alternative
#' @template distr
#' @param pValueFUN A function that calculates the p-value. Must have an argument \code{mu} that specifies the null-hypothesis. 
#' Defaults to \code{\link[hMean]{hMeanChiSqMu}}.
#' @template heterogeneity
#' @param wGamma Numeric vector of length \code{unique(thetahat) - 1} specifying weights used to
#' summarize the gamma values, i.e.,
#' the local minima of the p-value function between the thetahats. Default is a vector of 1s.
#' @template check_inputs
#' @return Returns a list containing confidence interval(s)
#' obtained by inverting the harmonic mean chi-squared test based on study-specific
#' estimates and standard errors. The list contains:
#' \item{CI}{Confidence interval(s).}\cr\cr
#' If the \code{alternative} is "none", the list also contains:
#' \item{gamma}{Local minima of the p-value function between the thetahats.}
#' \item{gammaMean}{Mean of all gammas weighted by \code{wGamma}.}
#' \item{gammaHMean}{Harmonic mean of all gammas weighted by \code{wGamma}.}
#' @export
#' @importFrom stats uniroot optimize qnorm weighted.mean
#' @importFrom methods formalArgs
hMeanChiSqCI <- function(thetahat, se, 
                         w = rep(1, length(thetahat)),
                         phi = NULL,
                         tau2 = NULL,
                         level = 0.95, 
                         alternative = "none",
                         distr = c("chisq", "f"),
                         pValueFUN = hMeanChiSqMu,
                         heterogeneity = c("additive", "multiplicative"),
                         wGamma = rep(1, length(unique(thetahat)) - 1),
                         check_inputs = TRUE){
  
  # Check inputs
  if(check_inputs){
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
              
              is.null(phi) || is.numeric(phi) & length(phi) == 1L & is.finite(phi) & 0 <= phi, # should phi be allowed to be < 1?
              is.null(tau2) || is.numeric(tau2) & length(tau2) == 1L & is.finite(tau2) & 0 <= tau2,
              !is.null(tau2) || !is.null(phi), # one of both must be given
              
              !is.null(alternative),
              length(alternative) == 1L,
              alternative %in% c("greater", "less", "two.sided", "none"),
              
              is.numeric(level),
              length(level) == 1L,
              is.finite(level),
              level > 0 & level < 1,
              
              is.numeric(wGamma),
              length(wGamma) == length(unique(thetahat)) - 1,
              
              !is.null(distr),
              length(distr) == 1L,
              
              !is.null(heterogeneity),
              length(heterogeneity) == 1L,
              
              (heterogeneity == "additive" & !is.null(tau2)) || (heterogeneity == "multiplicative" & !is.null(phi)))
  }
  
  
  # estimate heterogeneity
  distr <- match.arg(distr, several.ok = FALSE)
  if(length(se) == 1L) se <- rep(se, length(thetahat))
  
  # target function to compute the limits of the CI
  args <- alist(thetahat = thetahat, se = se, w = w, phi = phi, tau2 = tau2, mu = limit,
                alternative = alternative, distr = distr, heterogeneity = heterogeneity,
                bound = FALSE, check_inputs = FALSE)
  args <- args[names(args) %in% methods::formalArgs(pValueFUN)]
  target <- function(limit){
    do.call(`pValueFUN`, args) - alpha
  }
  # target <- function(limit){
  #   hMeanChiSqMu(thetahat = thetahat, se = se, w = w, mu = limit, phi = phi, tau2 = tau2,
  #                alternative = alternative, distr = distr, heterogeneity = heterogeneity,
  #                bound = FALSE, check_inputs = FALSE) - alpha
  # }
  
  ## sort 'thetahat', 'se', 'w'
  indOrd <- order(thetahat)
  thetahat <- thetahat[indOrd]; se <- se[indOrd]; w <- w[indOrd]
  
  ## minima are only searched between distinct thetahat elements
  thetahatUnique <- unique(thetahat)
  nThetahatUnique <- length(thetahatUnique)
  
  mini <- which.min(thetahat)
  maxi <- which.max(thetahat)
  mint <- thetahat[mini]
  maxt <- thetahat[maxi]
  minse <- se[mini]
  maxse <- se[maxi]
  alpha <- 1 - level
  z1 <- max(-stats::qnorm(alpha), 1)
  eps <- 1e-6
  factor <- 5
  if(alternative == "none"){
    
    ## ----------------------------
    ## find lower bound such that: lower < thetahat[1] AND target(lower) < 0 
    lower <- mint - z1 * minse
    while(target(lower) > 0){
      lower <- lower - minse
    }
    ## find root between 'lower' and 'thetahat[1]'
    CIlower <- stats::uniroot(f = target, lower = lower, upper = thetahat[1])$root
    
    ## -------------------------
    ## check between thetahats whether 'target' goes below 'alpha'
    ## if so, search CI limits
    CImiddle <- matrix(NA_real_, nrow = 2L, ncol = nThetahatUnique - 1L)
    gam <- matrix(NA_real_, nrow = nThetahatUnique - 1L, ncol = 2L)
    colnames(gam) <- c("minimum", "pvalue_fun/gamma")
    for(i in seq_len(nThetahatUnique - 1L)){
      opt <- optimize(f = target, lower = thetahatUnique[i],
                      upper = thetahatUnique[i + 1L])
      gam[i,] <- c(opt$minimum, opt$objective + alpha)
      if(opt$objective <= 0){
        CImiddle[1L, i] <- stats::uniroot(f = target, lower = thetahatUnique[i],
                                          upper = opt$minimum)$root
        CImiddle[2L, i] <- stats::uniroot(f = target, lower = opt$minimum,
                                          upper = thetahatUnique[i + 1L])$root
      }
    }
    CImiddle <- CImiddle[!is.na(CImiddle)]
    
    ## -------------------------
    ## find upper bound such that:
    ## upper > thetahat[length(thetahat)] AND target(upper) < 0 
    upper <- maxt + maxse
    while(target(upper) > 0){
      upper <- upper + z1 * maxse
    }
    ## find root between 'lower' and 'thetahat[1]'
    CIupper <- stats::uniroot(f = target, lower = thetahat[length(thetahat)],
                              upper = upper)$root
    CI <- matrix(c(CIlower, CImiddle, CIupper), ncol = 2, byrow = TRUE)
    colnames(CI) <- c("lower", "upper")
    return(list(CI = CI,
                gamma = gam,
                gammaMean = stats::weighted.mean(x = gam[,"pvalue_fun/gamma"], w = wGamma),
                gammaHMean = sum(wGamma) / sum(wGamma / gam[,"pvalue_fun/gamma"])))
    
  } else if(alternative == "two.sided"){
    lower <- stats::uniroot(f = target,
                            lower = mint - factor * z1 * minse,
                            upper = mint - eps * minse)$root
    upper <- stats::uniroot(f = target, lower = maxt + eps * maxse,
                            upper = maxt + factor * z1 * maxse)$root
    return(list(CI = cbind(lower, upper)))
    
  } else if(alternative == "greater"){
    lower <- stats::uniroot(f = target,
                            lower = mint - factor * z1 * minse,
                            upper = mint - eps * minse)$root
    upper <- Inf
    return(list(CI = cbind(lower, upper)))
    
  } else if(alternative == "less"){
    lower <- -Inf
    upper <- stats::uniroot(f = target,
                            lower = maxt + eps * maxse,
                            upper = maxt + factor * z1 * maxse)$root
    return(list(CI = cbind(lower, upper)))
    
  } else {
    stop(paste0("The argument 'alternative = ", alternative, "' has not been implemented yet."))
  }
  
}
