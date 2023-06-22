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
#' @template pValueFUN_args
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
    get_CI(
        alternative = alternative,
        f = f,
        bounds = bounds,
        wGamma = wGamma
    )
}

################################################################################
# Helper function that returns a function to optimize                          #
################################################################################

make_function <- function(
    thetahat,
    se,
    alpha,
    pValueFUN,
    pValueFUN_args
) {
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
    list(CI = cbind("lower" = -Inf, "upper" = upper))
}

get_CI_greater <- function(f, bounds) {
    lower <- stats::uniroot(
        f = f,
        lower = bounds[1],
        upper = bounds[2]
    )$root
    list(CI = cbind("lower" = lower, "upper" = Inf))
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
        list(CI = cbind("lower" = lower, "upper" = upper))
    )
}

get_CI_none <- function(f, bounds, wGamma) {

    # get the thetahats and save the lower and upper CI bounds
    idx <- c(1L, length(bounds))
    thetahat <- bounds[-idx]
    l <- bounds[idx[1L]]
    u <- bounds[idx[2L]]

    # Get the number of intervals between the thetahats
    n_intervals <- length(thetahat) - 1L

    # For the intervals in the middle, compute the minimum and the corresponding
    # p-value
    gam <- t(
        vapply(
            seq_len(n_intervals),
            function(i) {
                opt <- stats::optimize(
                    f = f,
                    lower = thetahat[i],
                    upper = thetahat[i + 1L]
                )
                c(opt$minimum, opt$objective)
            },
            double(2L)
        )
    )
    colnames(gam) <- c("minimum", "pvalue_fun/gamma")

    # Whereever the p-value function is negative at the minimum,
    # search for the two roots. Also add the lower and upper bound
    minima <- gam[, 2L]
    if (any(minima < 0)) {
        negative <- seq_len(nrow(gam))[minima < 0]
        min <- gam[, 1L][negative]
        CI <- vapply(
            seq_along(negative),
            function(i) {
                l <- stats::uniroot(
                    f = f,
                    lower = thetahat[negative[i]],
                    upper = min[i]
                )$root
                u <- stats::uniroot(
                    f = f,
                    lower = min[i],
                    upper = thetahat[negative[i] + 1L]
                )$root
                c(l, u)
            },
            double(2L)
        )
        CI <- matrix(c(l, CI, u), ncol = 2L, byrow = TRUE)
    } else {
        CI <- matrix(c(l, u), ncol = 2L, byrow = TRUE)
    }
    colnames(CI) <- c("lower", "upper")

    # Increase the y-coordinate of the minima by alpha
    gam[, 2L] <- gam[, 2L] + get("alpha", envir = environment(f))

    # return
    list(
        CI = CI,
        gamma = gam,
        gammaMean = stats::weighted.mean(
            x = gam[, 2L],
            w = wGamma
        ),
        gammaHMean = sum(wGamma) / sum(wGamma / gam[, 2L])
    )
}

get_CI <- function(alternative, f, bounds, wGamma) {
    switch(
        alternative,
        "none" = get_CI_none(f = f, bounds = bounds, wGamma),
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
    while (f(lower) > 0) {
        lower <- lower - minse
    }

    ## find the root
    lower <- stats::uniroot(f, lower = lower, upper = mint)$root

    ## find upper bound such that:
    ## upper > thetahat[length(thetahat)] AND target(upper) < 0
    upper <- maxt + maxse
    while (f(upper) > 0) {
        upper <- upper + z1 * maxse
    }

    ## find the root
    upper <- stats::uniroot(f = f, lower = maxt, upper = upper)$root

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
