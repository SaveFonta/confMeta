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
#' @examples
#' n <- 15
#' mean <- 0
#' sd <- 1.1
#' shape <- 5
#' rate <- 5
#' thetahat <- rnorm(n, mean = mean, sd = sd)
#' se <- rgamma(n, shape = shape, rate = rate)
#' heterogeneity <- "none"
#' phi <- if (heterogeneity == "multiplicative") {
#'     estimatePhi(thetahat = thetahat, se = se)
#' } else {
#'     NULL
#' }
#' tau2 <- if (heterogeneity == "additive") {
#'     estimateTau2(thetahat = thetahat, se = se)
#' } else {
#'     NULL
#' }
#' mu <- seq(
#'   min(thetahat) - 0.5 * max(se),
#'   max(thetahat) + 0.5 * max(se),
#'   length.out = 1e5
#' )
#' alpha <- 0.05
#' funs <- list(
#'     "pearson" = hMean::pPearsonMu,
#'     "hMean" = hMean::hMeanChiSqMu,
#'     "k-Trials" = hMean::kTRMu,
#'     "edgington" = hMean::pEdgingtonMu,
#'     "fisher" = hMean::pFisherMu
#' )
#' p_vals <- do.call(
#'     "cbind",
#'     lapply(
#'         funs,
#'         function(f) {
#'             f(
#'                 thetahat = thetahat,
#'                 mu = mu,
#'                 se = se,
#'                 heterogeneity = heterogeneity,
#'                 phi = phi,
#'                 tau2 = tau2,
#'                 check_inputs = FALSE
#'             )
#'         }
#'     )
#' )
#' cis <- lapply(
#'     funs,
#'     function(f) {
#'         hMeanChiSqCI(
#'             thetahat = thetahat,
#'             se = se,
#'             level = 1 - alpha,
#'             alternative = "none",
#'             pValueFUN = f,
#'             pValueFUN_args = list(
#'                 check_inputs = FALSE,
#'                 heterogeneity = heterogeneity,
#'                 phi = phi,
#'                 tau2 = tau2
#'             )
#'         )
#'     }
#' )
#' plot_res <- function(
#'     mu,
#'     p_vals,
#'     cis = NULL,
#'     barheight = 0.05
#' ) {
#'     opar <- par(no.readonly = TRUE)
#'     par(las = 1)
#'     matplot(
#'         mu,
#'         p_vals,
#'         type = "l", lty = 1, lwd = 3,
#'         ylab = "p-value function", xlab = expression(mu)
#'     )
#'     legend("topleft",
#'         col = c(1, 2, 3, 4, 5),
#'         lwd = 3,
#'         lty = 1,
#'         legend = c("Pearson", "hMean", "k-Trials", "Edgington", "Fisher"),
#'         bty = "n",
#'         cex = 2
#'     )
#'     abline(h = 0.05, lty = 2)
#'     if (!is.null(cis)) {
#'         cis <- lapply(cis, "[[", i = "CI")
#'         jitter_inc <- 0.001
#'         jitter <- jitter_inc
#'         for (j in seq_along(cis)) {
#'             x <- cis[[j]]
#'             y_horiz <- alpha + (-1)^j * jitter
#'             if (j %% 2 == 0L) jitter <- jitter + jitter_inc
#'             for (i in seq_len(nrow(x))) {
#'                 l <- x[i, "lower"]
#'                 u <- x[i, "upper"]
#'                 lty <- 1
#'                 lwd <- 2
#'                 segments( # horizontal
#'                     x0 = l,
#'                     x1 = u,
#'                     y0 = y_horiz,
#'                     y1 = y_horiz,
#'                     lty = lty,
#'                     lwd = lwd,
#'                     col = j
#'                 )
#'                 segments( # error bar left
#'                     x0 = l,
#'                     x1 = l,
#'                     y0 = y_horiz - barheight / 2,
#'                     y1 = y_horiz + barheight / 2,
#'                     lty = lty,
#'                     lwd = lwd,
#'                     col = j
#'                 )
#'                 segments( # error bar left
#'                     x0 = u,
#'                     x1 = u,
#'                     y0 = y_horiz - barheight / 2,
#'                     y1 = y_horiz + barheight / 2,
#'                     lty = lty,
#'                     lwd = lwd,
#'                     col = j
#'                 )
#'             }
#'         }
#'     }
#'     par(opar)
#' }
#' plot_res(mu = mu, p_vals = p_vals, cis = cis)
#'
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
    # This is calls the p-value function with specified
    # args and subtracts alpha
    f <- make_function(
        thetahat = thetahat,
        se = se,
        alpha = 1 - level,
        pValueFUN = pValueFUN,
        pValueFUN_args = pValueFUN_args
    )

    # sort thetahat and se
    o <- order(thetahat, decreasing = FALSE)
    thetahat <- thetahat[o]
    se <- se[o]

    # Check if CI even exists: This is the case if
    # the function f(thetahat) returns at least one
    # positive value
    f_thetahat <- f(thetahat)
    if (all(f_thetahat <= 0)) {
        # If it does not exist, return same format but all NAs
        out <- list(
            CI = matrix(rep(NA_real_, 2L), ncol = 2L),
            gamma = matrix(rep(NA_real_, 2L), ncol = 2L),
            gammaMean = NA_real_,
            gammaHMean = NA_real_
        )
        colnames(out$CI) <- c("lower", "upper")
        colnames(out$gamma) <- c("minimum", "pvalue_fun/gamma")
    } else {
        # If the CI does exist:
        # 1. Determine the smallest and largest thetahat where f(thetahat) > 0
        # 2. Find the lower and upper bounds based on these thetahats
        # 3. Corners/Cusps are always at thetahats. Thus, we search between
        #    the lower bound, thetahat_min, all the thetahats in between and
        #    finally thetahat_max and the upper bound, this is implemented in
        #    the function get_CI

        # 1.
        thetahat_pos <- which(f_thetahat > 0)
        idx_min <- min(thetahat_pos)
        idx_max <- max(thetahat_pos)
        thetahat_min <- thetahat[idx_min]
        se_min <- se[idx_min]
        thetahat_max <- thetahat[idx_max]
        se_max <- se[idx_max]
        # 2.
        lower <- find_lower(
            f = f,
            thetahat_min = thetahat_min,
            se_min = se_min
        )
        upper <- find_upper(
            f = f,
            thetahat_max = thetahat_max,
            se_max = se_max
        )
        # 3.
        ## Get the thetahats we need to examine. These are all in
        ## between thetahat_min and thetahat_max
        thetahat <- thetahat[idx_min:idx_max]

        ## Get the number of intervals between these thetahats
        n_intervals <- length(thetahat) - 1L

        ## For the intervals in the middle, compute the minimum and the
        ## corresponding p-value
        if (n_intervals == 0) {
            gam <- matrix(NA_real_, ncol = 2L, nrow = 1L)
        } else {
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
        }
        colnames(gam) <- c("minimum", "pvalue_fun/gamma")

        # Whereever the p-value function is negative at the minimum,
        # search for the two roots. Also add the lower and upper bound
        # If there is no minimum (i.e. only one thetahat is positive),
        # then, we can also just use lower & upper for the CI
        minima <- gam[, 2L]
        one_pos_theta_only <- length(minima) == 1L && is.na(minima)
        exist_neg_minima <- any(minima < 0)
        search_roots <- !one_pos_theta_only && exist_neg_minima
        if (search_roots) {
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
            CI <- matrix(c(lower, CI, upper), ncol = 2L, byrow = TRUE)
        } else {
            CI <- matrix(c(lower, upper), ncol = 2L, byrow = TRUE)
        }
        colnames(CI) <- c("lower", "upper")

        # Increase the y-coordinate of the minima by alpha
        if (!one_pos_theta_only) {
            gam[, 2L] <- gam[, 2L] + get("alpha", envir = environment(f))
        }

        # return
        out <- list(
            CI = CI,
            gamma = gam,
            gammaMean = mean(gam[, 2L]),
            # gammaMean = stats::weighted.mean(
            #     x = gam[, 2L],
            #     w = wGamma
            # ),
            # gammaHMean = sum(wGamma) / sum(wGamma / gam[, 2L])
            gammaHMean = nrow(gam) / sum(nrow(gam) / gam[, 2L])
        )
    }
    out
}


################################################################################
# Helper function that returns a function to optimize                          #
################################################################################

find_lower <- function(thetahat_min, se_min, f) {
    lower <- thetahat_min - se_min
    while (f(lower) > 0) {
        lower <- lower - se_min
    }
    stats::uniroot(
        f = f,
        lower = lower,
        upper = thetahat_min
    )$root
}

find_upper <- function(thetahat_max, se_max, f) {
    upper <- thetahat_max + se_max
    while (f(upper) > 0) {
        upper <- upper + se_max
    }
    stats::uniroot(
        f = f,
        lower = thetahat_max,
        upper = upper
    )$root
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
#
# get_CI_less <- function(f, bounds) {
#     upper <- stats::uniroot(
#         f = f,
#         lower = bounds[1],
#         upper = bounds[2]
#     )$root
#     list(CI = cbind("lower" = -Inf, "upper" = upper))
# }
#
# get_CI_greater <- function(f, bounds) {
#     lower <- stats::uniroot(
#         f = f,
#         lower = bounds[1],
#         upper = bounds[2]
#     )$root
#     list(CI = cbind("lower" = lower, "upper" = Inf))
# }
#
# get_CI_twosided <- function(f, bounds) {
#     lower <- stats::uniroot(
#         f = f,
#         lower = bounds[1],
#         upper = bounds[2]
#     )$root
#     upper <- stats::uniroot(
#         f = f,
#         lower = bounds[3],
#         upper = bounds[4]
#     )$root
#     return(
#         list(CI = cbind("lower" = lower, "upper" = upper))
#     )
# }
#
# get_CI_none <- function(f, bounds, wGamma) {
#
#     # get the thetahats and save the lower and upper CI bounds
#     idx <- c(1L, length(bounds))
#     thetahat <- bounds[-idx]
#     l <- bounds[idx[1L]]
#     u <- bounds[idx[2L]]
#
#     # Get the number of intervals between the thetahats
#     n_intervals <- length(thetahat) - 1L
#
#     # For the intervals in the middle, compute the minimum and the
#     # corresponding p-value
#     gam <- t(
#         vapply(
#             seq_len(n_intervals),
#             function(i) {
#                 opt <- stats::optimize(
#                     f = f,
#                     lower = thetahat[i],
#                     upper = thetahat[i + 1L]
#                 )
#                 c(opt$minimum, opt$objective)
#             },
#             double(2L)
#         )
#     )
#     colnames(gam) <- c("minimum", "pvalue_fun/gamma")
#
#     # Whereever the p-value function is negative at the minimum,
#     # search for the two roots. Also add the lower and upper bound
#     minima <- gam[, 2L]
#     if (any(minima < 0)) {
#         negative <- seq_len(nrow(gam))[minima < 0]
#         min <- gam[, 1L][negative]
#         CI <- vapply(
#             seq_along(negative),
#             function(i) {
#                 l <- stats::uniroot(
#                     f = f,
#                     lower = thetahat[negative[i]],
#                     upper = min[i]
#                 )$root
#                 u <- stats::uniroot(
#                     f = f,
#                     lower = min[i],
#                     upper = thetahat[negative[i] + 1L]
#                 )$root
#                 c(l, u)
#             },
#             double(2L)
#         )
#         CI <- matrix(c(l, CI, u), ncol = 2L, byrow = TRUE)
#     } else {
#         CI <- matrix(c(l, u), ncol = 2L, byrow = TRUE)
#     }
#     colnames(CI) <- c("lower", "upper")
#
#     # Increase the y-coordinate of the minima by alpha
#     gam[, 2L] <- gam[, 2L] + get("alpha", envir = environment(f))
#
#     # return
#     list(
#         CI = CI,
#         gamma = gam,
#         gammaMean = stats::weighted.mean(
#             x = gam[, 2L],
#             w = wGamma
#         ),
#         gammaHMean = sum(wGamma) / sum(wGamma / gam[, 2L])
#     )
# }
#
# get_CI <- function(alternative, f, bounds, wGamma) {
#     switch(
#         alternative,
#         "none" = get_CI_none(f = f, bounds = bounds, wGamma),
#         "less" = get_CI_less(f = f, bounds = bounds),
#         "greater" = get_CI_greater(f = f, bounds = bounds),
#         "two.sided" = get_CI_twosided(f = f, bounds = bounds)
#     )
# }
#
################################################################################
# Helper functions to calculate the bounds                                     #
################################################################################
#
# get_bounds_none <- function(f, thetahat, mint, maxt, minse, maxse, z1) {
#
#     ## find lower bound such that: lower < thetahat[1] AND target(lower) < 0
#     lower <- mint - z1 * minse
#     while (f(lower) > 0) {
#         lower <- lower - minse
#     }
#
#     ## find the root
#     lower <- stats::uniroot(f, lower = lower, upper = mint)$root
#
#     ## find upper bound such that:
#     ## upper > thetahat[length(thetahat)] AND target(upper) < 0
#     upper <- maxt + maxse
#     while (f(upper) > 0) {
#         upper <- upper + z1 * maxse
#     }
#
#     ## find the root
#     upper <- stats::uniroot(f = f, lower = maxt, upper = upper)$root
#
#     c(lower, upper)
# }
#
# get_bounds <- function(alternative, thetahat, se, level, f) {
#
#     ## sort 'thetahat', 'se'
#     o <- order(thetahat)
#     thetahat <- thetahat[o]
#     se <- se[o]
#
#     ## minima are only searched between distinct thetahat elements
#     thetahat <- unique(thetahat)
#
#     # get the lower and upper bound
#     mini <- which.min(thetahat)
#     maxi <- which.max(thetahat)
#     mint <- thetahat[mini]
#     maxt <- thetahat[maxi]
#     minse <- se[mini]
#     maxse <- se[maxi]
#
#     # set some constants
#     alpha <- 1 - level
#     z1 <- max(-stats::qnorm(alpha), 1)
#     # eps <- 1e-6
#     # factor <- 5
#
#     if (alternative == "none") {
#         b <- get_bounds_none(
#             f = f,
#             mint = mint,
#             minse = minse,
#             maxt = maxt,
#             maxse = maxse,
#             z1 = z1
#         )
#     }
#
#     # return bounds
#     switch(
#         alternative,
#         "none" = c(
#             b[1], thetahat, b[2]
#         )#,
#         # "less" = c(
#         #     maxt + eps * maxse,
#         #     maxt + factor * z1 * maxse
#         # ),
#         # "greater" = c(
#         #     mint - factor * z1 * minse,
#         #     mint - eps * minse
#         # ),
#         # "two.sided" = c(
#         #     mint - factor * z1 * minse,
#         #     mint - eps * minse,
#         #     maxt + eps * maxse,
#         #     maxt + factor * z1 * maxse
#         # )
#     )
# }
