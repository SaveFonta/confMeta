#' Calculate confidence intervals based on the harmonic mean chi-squared test
#' 
#' @template thetahat
#' @template se
#' @param level Numeric vector of length 1 specifying the level of the
#' confidence interval. Defaults to 0.95.
#' @template alternative
#' @param pValueFUN A function that calculates the p-value. Must have arguments
#' \code{thetahat} and \code{se} as these are passed by this function.
#' Must further have an argument \code{mu} that specifies the null-hypothesis.
#' Defaults to \code{\link[hMean]{hMeanChiSqMu}}.
#' @param wGamma Numeric vector of length \code{unique(thetahat) - 1} specifying
#' weights used to
#' summarize the gamma values, i.e., the local minima of the p-value function
#' between the thetahats. Default is a vector of 1s.
#' @template check_inputs
#' @param pValueFUN_args A named \code{list} with arguments passed to
#' \code{pValueFUN}. This list must contains all arguments of \code{pValueFUN}
#' that do not have a Default. Arguments \code{thetahat}, \code{se}, and
#' \code{mu} are automatically passed to \code{pValueFUN} and must thus not be
#' included in this named list.
#' @return Returns a list containing confidence interval(s)
#' obtained by inverting the harmonic mean chi-squared test based on
#' study-specific estimates and standard errors. The list contains:
#' \item{CI}{Confidence interval(s).}\cr\cr
#' If the \code{alternative} is "none", the list also contains:
#' \item{gamma}{Local minima of the p-value function between the thetahats.}
#' \item{gammaMean}{Mean of all gammas weighted by \code{wGamma}.}
#' \item{gammaHMean}{Harmonic mean of all gammas weighted by \code{wGamma}.}
#' @export
hMeanChiSqCI <- function(
  thetahat,
  se,
  level = 0.95, 
  alternative = "none",
  wGamma = rep(1, length(unique(thetahat)) - 1),
  check_inputs = TRUE,
  pValueFUN = hMeanChiSqMu,
  pValueFUN_args
) {
  
  # Check inputs
  if (check_inputs) {
    check_inputs_CI(
      thetahat = thetahat,
      se = se,
      level = level,
      alternative = alternative,
      wGamma = wGamma,
      check_inputs = check_inputs,
      pValueFUN = pValueFUN,
      pValueFUN_args = pValueFUN_args
    )
  }
  
  # expand se
  if (length(se) == 1L) se <- rep(se, length(thetahat))

  # target function to compute the limits of the CI
  ## pass dotargs, thetahat, se
  args <- append(pValueFUN_args, list(thetahat = thetahat, se = se))
  ## add mu
  args <- append(args, alist(mu = limit))
  ## For the remaining arguments, use the defaults
  args <- append(
    args,
    formals(pValueFUN)[!methods::formalArgs(pValueFUN) %in% names(args)]
  )
  ## define target function
  target <- function(limit) {
    do.call(`pValueFUN`, args) - alpha
  }
  
  ## sort 'thetahat', 'se', 'w'
  indOrd <- order(thetahat)
  thetahat <- thetahat[indOrd]
  se <- se[indOrd]
  
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
  if (alternative == "none") {

    ## ----------------------------
    ## find lower bound such that: lower < thetahat[1] AND target(lower) < 0 
    lower <- mint - z1 * minse
    while(target(lower) > 0) {
      lower <- lower - minse
    }
    ## find root between 'lower' and 'thetahat[1]'
    CIlower <- stats::uniroot(
      f = target,
      lower = lower,
      upper = thetahat[1]
    )$root
    
    ## -------------------------
    ## check between thetahats whether 'target' goes below 'alpha'
    ## if so, search CI limits
    CImiddle <- matrix(NA_real_, nrow = 2L, ncol = nThetahatUnique - 1L)
    gam <- matrix(NA_real_, nrow = nThetahatUnique - 1L, ncol = 2L)
    colnames(gam) <- c("minimum", "pvalue_fun/gamma")
    for(i in seq_len(nThetahatUnique - 1L)) {
      opt <- stats::optimize(
        f = target,
        lower = thetahatUnique[i],
        upper = thetahatUnique[i + 1L]
      )
      gam[i,] <- c(opt$minimum, opt$objective + alpha)
      if (opt$objective <= 0) {
        CImiddle[1L, i] <- stats::uniroot(
          f = target,
          lower = thetahatUnique[i],
          upper = opt$minimum
        )$root
        CImiddle[2L, i] <- stats::uniroot(
          f = target,
          lower = opt$minimum,
          upper = thetahatUnique[i + 1L]
        )$root
      }
    }
    CImiddle <- CImiddle[!is.na(CImiddle)]
    
    ## -------------------------
    ## find upper bound such that:
    ## upper > thetahat[length(thetahat)] AND target(upper) < 0 
    upper <- maxt + maxse
    while (target(upper) > 0) {
      upper <- upper + z1 * maxse
    }
    ## find root between 'lower' and 'thetahat[1]'
    CIupper <- stats::uniroot(
      f = target,
      lower = thetahat[length(thetahat)],
      upper = upper
    )$root
    CI <- matrix(c(CIlower, CImiddle, CIupper), ncol = 2, byrow = TRUE)
    colnames(CI) <- c("lower", "upper")
    return(
      list(
        CI = CI,
        gamma = gam,
        gammaMean = stats::weighted.mean(
          x = gam[,"pvalue_fun/gamma"],
          w = wGamma
        ),
        gammaHMean = sum(wGamma) / sum(wGamma / gam[,"pvalue_fun/gamma"])
      )
    )
    
  } else if (alternative == "two.sided") {
    lower <- stats::uniroot(
      f = target,
      lower = mint - factor * z1 * minse,
      upper = mint - eps * minse
    )$root
    upper <- stats::uniroot(
      f = target,
      lower = maxt + eps * maxse,
      upper = maxt + factor * z1 * maxse
    )$root
    return(
      list(CI = cbind(lower, upper))
    )
    
  } else if (alternative == "greater") {
    lower <- stats::uniroot(
      f = target,
      lower = mint - factor * z1 * minse,
      upper = mint - eps * minse
    )$root
    upper <- Inf
    return(list(CI = cbind(lower, upper)))
    
  } else if (alternative == "less") {
    lower <- -Inf
    upper <- stats::uniroot(
      f = target,
      lower = maxt + eps * maxse,
      upper = maxt + factor * z1 * maxse
    )$root
    return(list(CI = cbind(lower, upper)))
    
  } else {
    stop(
      paste0(
        "The argument 'alternative = ",
        alternative,
        "' has not been implemented yet."
      )
    )
  }
}


################################################################################
hMeanChiSqCI_new <- function(
  thetahat,
  se,
  level = 0.95, 
  alternative = "none",
  wGamma = rep(1, length(unique(thetahat)) - 1),
  check_inputs = TRUE,
  pValueFUN = hMeanChiSqMu,
  pValueFUN_args
) {
  
  if (check_inputs) {
    check_inputs_CI(
      thetahat = thetahat,
      se = se,
      level = level,
      alternative = alternative,
      wGamma = wGamma,
      pValueFUN = pValueFUN,
      pValueFUN_args = pValueFUN_args
    )
  }
  
  # Get the function we need to optimise
  f <- make_function(
    thetahat = thetahat,
    se = se,
    alpha = 1 - level,
    pValueFUN = pValueFUN,
    pValueFUN_args = pValueFUN_args
  )
  
  # Get the bounds
  bounds <- get_bounds(
    alternative = alternative,
    thetahat = thetahat,
    se = se,
    level = level,
    f = f
  )
  
  # Find the CIs
  CIs <- get_CI(
    
  )
  
    
  
  
}

################################################################################
# Helper function that returns a function to optimize                          #
################################################################################

make_function <- function(thetahat, se, alpha, pValueFUN, pValueFUN_args){
  # Add/Overwrite thetahat and se args
  pValueFUN_args$thetahat <- thetahat
  pValueFUN_args$se <- se
  # Add mu argument
  if ("mu" %in% names(pValueFUN_args)) pValueFUN_args$mu <- NULL
  pValueFUN_args <- append(pValueFUN_args, alist(mu = limit))
  ## For the remaining arguments, use the defaults
  forms <- formals(pValueFUN)
  nforms <- names(formals)
  pValueFUN_args <- append(
    pValueFUN_args,
    forms[!nforms %in% names(pValueFUN_args)]
  )
  ## Check whether all arguments are there
  available_args <- nforms %in% names(pValueFUN_args)
  if (!all(available_args)) {
    stop(
      paste0(
        "List pValueFUN_args is missing argument(s) '", 
        paste0(nforms[!available_args], collapse = "', '"),
        "'."
      )
    )
  }
  
  function(limit) {
    do.call(pValueFUN, pValueFUN_args) - alpha
  }
}

################################################################################
# Helper functions to calculate the CIs                                        #
################################################################################

get_CI_less <- function(f, bounds) {
    upper <- stats::uniroot(
      f = f,
      lower = bounds[1],
      upper = bounds[2]
    )$root
    list(CI = cbind(-Inf, upper))
}

get_CI_greater <- function(f, bounds) {
    lower <- stats::uniroot(
      f = f,
      lower = bounds[1],
      upper = bounds[2]
    )$root
    list(CI = cbind(lower, -Inf))
}

get_CI_twosided <- function(f, bounds) {
    lower <- stats::uniroot(
      f = f,
      lower = bounds[1],
      upper = bounds[2]
    )$root
    upper <- stats::uniroot(
      f = f,
      lower = bounds[3],
      upper = bounds[4]
    )$root
    return(
      list(CI = cbind(lower, upper))
    )
}

get_CI_none <- function(f, bounds) {
    n_intervals <- length(bounds) - 1L
    CI <- gam <- matrix(NA_real_, ncol = 2L, nrow = n_intervals)
    colnames(gam) <- c("minimum", "pvalue_fun/gamma")
    colnames(CI) <- c("lower", "upper")
    for (i in seq_len(n_intervals)) {
      if (i == 1L || i == n_intervals) {
        get_root_ends <- function(f, bounds) {
          
        }
      } else {
        opt <- stats::optimize(
          f = f,
          lower = bounds[i],
          upper = bounds[i + 1L]
        )
        mini <- opt$minimum
        obj <- opt$objective
        gam[i, ] <- c(mini, obj + alpha)
        if (obj <= 0) {
          CI[i, 1L] <- stats::uniroot(
            f = f,
            lower = bounds[i],
            upper = mini
          )$root
          CI[i, 2L] <- stats::uniroot(
            f = target,
            lower = mini
            upper = bounds[i + 1L]
          )$root
        }
      }
    }
}

get_CI <- function(alternative, f, bounds, alpha) {
  switch(
    alternative,
    "none" = get_CI_none(f = f, bounds = bounds, alpha = alpha),
    "less" = get_CI_less(f = f, bounds = bounds),
    "greater" = get_CI_greater(f = f, bounds = bounds),
    "two.sided" = get_CI_twosided(f = f, bounds = bounds)
  )
}

################################################################################
# Helper functions to calculate the bounds                                     #
################################################################################

get_bounds_none <- function(f, thetahat, mint, maxt, minse, maxse, z1) {
  
    ## find lower bound such that: lower < thetahat[1] AND target(lower) < 0 
    lower <- mint - z1 * minse
    while(f(lower) > 0) {
      lower <- lower - minse
    }
    
    ## find upper bound such that:
    ## upper > thetahat[length(thetahat)] AND target(upper) < 0 
    upper <- maxt + maxse
    while (f(upper) > 0) {
      upper <- upper + z1 * maxse
    }
    
    c(lower, upper)
}

get_bounds <- function(alternative, thetahat, se, level, f) {
  
    ## sort 'thetahat', 'se'
    o <- order(thetahat)
    thetahat <- thetahat[o]
    se <- se[o]
    
    ## minima are only searched between distinct thetahat elements
    thetahat <- unique(thetahat)
    
    # get the lower and upper bound
    mini <- which.min(thetahat)
    maxi <- which.max(thetahat)
    mint <- thetahat[mini]
    maxt <- thetahat[maxi]
    minse <- se[mini]
    maxse <- se[maxi]
    
    # set some constants
    alpha <- 1 - level
    z1 <- max(-stats::qnorm(alpha), 1)
    eps <- 1e-6
    factor <- 5
    
    if (alternative == "none") {
      b <- get_bounds_none(
          f = f,
          mint = mint,
          minse = minse,
          maxt = maxt,
          maxse = maxse,
          z1 = z1
      )
    }
    
    # return bounds
    switch(
      alternative,
      "none" = c(
        b[1], thetahat, b[2]
      ),
      "less" = c(
        maxt + eps * maxse,
        maxt + factor * z1 * maxse
      ),
      "greater" = c(
        mint - factor * z1 * minse,
        mint - eps * minse
      ),
      "two.sided" = c(
        mint - factor * z1 * minse,
        mint - eps * minse,
        maxt + eps * maxse,
        maxt + factor * z1 * maxse
      )
    )
}
