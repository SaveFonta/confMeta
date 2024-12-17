get_ci <- function(
    estimates,
    SEs,
    conf_level,
    p_fun
) {

    # Keep a copy of estimates for later on since `estimates` is overwritten
    orig_est <- estimates
    orig_se <- SEs

    # Get the function we need to optimise
    # This is calls the p-value function with specified
    # args and subtracts alpha
    alpha <- 1 - conf_level
    f <- make_function(
        estimates = estimates,
        SEs = SEs,
        alpha = alpha,
        p_fun = p_fun
    )

    # remove duplicates and sort estimates and SEs
    keep <- !duplicated(estimates)
    estimates <- estimates[keep]
    SEs <- SEs[keep]
    o <- order(estimates, decreasing = FALSE)
    estimates <- estimates[o]
    SEs <- SEs[o]

    # Check if CI exists: This is the case if
    # the function f(estimates) returns at least one
    # positive value or we can find a local maximum x
    # between the estimates where f(x) > 0.
    # Also, keep track of the status:
    # - 0 = estimate
    # - 1 = maximum
    # - 2 = minimum
    estimates <- matrix(
        c(estimates, f(estimates), rep(0, length(estimates))),
        ncol = 3L,
        dimnames = list(NULL, c("x", "y", "status"))
    )
    ## search for local maxima in between estimates
    maxima <- find_optima(estimates = estimates[, 1L], f = f, maximum = TRUE)
    ## Find out which of these maxima is relevant, i.e. it has a higher p-value
    ## than both, the next smaller and the next larger estimate
    isRelevant_max <- is_relevant(
        f_estimates = estimates[, 2L],
        f_extremum = maxima[, 2L],
        maximum = TRUE
    )

    ## For searching CIs: Add the relevant maxima to the estimates
    ## Here, we only care about the maxima since they might be > 0
    ## even though none of the estimates are
    if (any(isRelevant_max)) {
        estimates <- rbind(estimates, maxima[isRelevant_max, ])
        ## Sort this by the x-coordinate
        o <- order(estimates[, 1L], decreasing = FALSE)
        estimates <- estimates[o, ]
    }
    f_estimates <- estimates[, 2L]
    estimates <- estimates[, 1L]

    # Calculate p_max
    idx <- f_estimates == max(f_estimates)
    p_max <- matrix(
        c(x = estimates[idx], y = f_estimates[idx] + alpha),
        ncol = 2L,
        dimnames = list(NULL, c("x", "y"))
    )

    # Calculate the AUCC and the ratio
    # Basically we need to integrate p_fun and for the ratio we need to
    # also integrate from -Infinity up to the x coordinate of p_max
    if (nrow(p_max) > 1L) {
        msg <- paste0(
            "More than one maximum of the p-value function found. ",
            "The AUCC ratio is thus "
        )
        warning(msg)
        aucc <- NA_real_
        aucc_ratio <- NA_real_
    } else {
        a <- integrate(
            p_fun,
            estimates = orig_est,
            SEs = orig_se,
            lower = -Inf,
            upper = p_max[, "x"]
        )$value
        b <- integrate(
            p_fun,
            estimates = orig_est,
            SEs = orig_se,
            lower = p_max[, "x"],
            upper = Inf
        )$value
        aucc <- a + b
        aucc_ratio <- if (a == 0) {
            1
        } else if (b == 0) {
            -1
        } else if (a == b) {
            0
        } else {
            (b - a) / (aucc)
        }
    }

    if (all(f_estimates <= 0)) {
        # If it does not exist, return same format
        out <- list(
            CI = matrix(rep(NA_real_, 2L), ncol = 2L),
            gamma = matrix(rep(NA_real_, 2L), ncol = 2L),
            p_max = p_max,
            p_0 = matrix(
                c(0, f(0) + alpha),
                ncol = 2L,
                dimnames = list(NULL, c("x", "y"))
            ),
            aucc = aucc,
            aucc_ratio = aucc_ratio
            # forest_plot_thetahat = thetahat[idx],
            # forest_plot_f_thetahat = f_thetahat[idx] + alpha
        )
        colnames(out$CI) <- c("lower", "upper")
        colnames(out$gamma) <- c("x", "y")
    } else {
        # If the CI does exist:
        # 1. Determine the smallest and largest estimate where f(estimate) > 0
        # 2. Find the lower and upper bounds based on these estimates
        # 3. Corners/Cusps are always at estimates. Thus, we search between
        #    the lower bound, estimate_min, all the estimates in between and
        #    finally estimate_max and the upper bound, this is implemented in
        #    the function get_CI

        # 1.
        estimates_pos <- which(f_estimates > 0)
        idx_min <- min(estimates_pos)
        idx_max <- max(estimates_pos)
        estimates_min <- estimates[idx_min]
        estimates_max <- estimates[idx_max]
        step <- max(SEs)

        # 2.
        lower <- find_lower(
            f = f,
            estimates_min = estimates_min,
            SEs_min = step
        )
        upper <- find_upper(
            f = f,
            estimates_max = estimates_max,
            SEs_max = step
        )
        # 3.
        ## Get the estimates we need to examine. These are all in
        ## between estimates_min and estimates_max.
        estimates <- estimates[idx_min:idx_max]
        f_estimates <- f_estimates[idx_min:idx_max]

        ## Get the number of intervals between these estimates
        n_intervals <- length(estimates) - 1L

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
                            lower = estimates[i],
                            upper = estimates[i + 1L]
                        )
                        c(opt$minimum, opt$objective)
                    },
                    double(2L)
                )
            )
        }
        colnames(gam) <- c("x", "y")

        # Whereever the p-value function is negative at the minimum,
        # search for the two roots. Also add the lower and upper bound
        # If there is no minimum (i.e. only one estimate is positive),
        # then, we can also just use lower & upper for the CI
        minima <- gam[, 2L]
        # only search roots if there is more than one positive f(estimate)
        # and there is at least one negative minimum
        one_pos_theta_only <- length(minima) == 1L && is.na(minima)
        exist_neg_minima <- any(minima < 0)
        search_roots <- !one_pos_theta_only && exist_neg_minima
        if (search_roots) {
            # Now that we know all the minima between the smallest and largest
            # positive estimate, we need to apply an algorithm to find the
            # roots. These exist if the minimum between two estimates i
            # and j is negative and both f(estimate[i]) and f(estimate[j])
            # are positive.

            # In order to find the correct intervals to search, we first
            # need to find all negative minima and then for each of them
            # the closest smaller estimate where f_estimate > 0 and the
            # closest larger estimate where f_estimate > 0. This is
            # is implemented in the functions get_search_interval and
            # find_closest_thetas
            intervals <- get_search_interval(
                x_max = estimates,
                y_max = f_estimates,
                x_min = gam[, 1L],
                y_min = gam[, 2L]
            )
            CI <- vapply(
                seq_len(ncol(intervals)),
                function(i) {
                    l <- stats::uniroot(
                        f = f,
                        lower = intervals[1L, i],
                        upper = intervals[3L, i]
                    )$root
                    u <- stats::uniroot(
                        f = f,
                        lower = intervals[3L, i],
                        upper = intervals[2L, i]
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
        # -> the if clause is only there because the object `gam`
        # does not exist if the p-value function is positive for
        # only one estimate/maximum
        if (!one_pos_theta_only) {
            gam[, 2L] <- gam[, 2L] + alpha
        }


        # return
        # Calculate p_max
        out <- list(
            CI = CI,
            gamma = gam,
            p_max = p_max,
            p_0 = matrix(
                c(0, f(0) + alpha),
                ncol = 2L,
                dimnames = list(NULL, c("x", "y"))
            ),
            aucc = aucc,
            aucc_ratio = aucc_ratio
            # gammaMean = mean(gam[, 2L]),
            # gammaHMean = nrow(gam) / sum(nrow(gam) / gam[, 2L]),
            # forest_plot_estimates = estimates,
            # forest_plot_f_estimates = f_estimates + alpha
        )
    }
    out
}

################################################################################
# Helper function to find out whether a local maximum is relevant or not       #
# Relevant here means that it has a higher p-value than the next smaller and   #
# the next larger effect estimate                                              #
################################################################################

is_relevant <- function(f_estimates, f_extremum, maximum) {
    lower <- f_estimates[-length(f_estimates)]
    upper <- f_estimates[-1L]
    if (maximum) {
        f_extremum > lower & f_extremum > upper
    } else {
        f_extremum < lower & f_extremum < upper
    }
}

################################################################################
# Helper function to find local maxima/minima between estimates                #
################################################################################

find_optima <- function(f, estimates, maximum, ...) {

    n_intervals <- length(estimates) - 1L
    status <- if (maximum) 1 else 2
    out <- t(
        vapply(
            seq_len(n_intervals),
            function(i) {
                opt <- stats::optimize(
                    f = f,
                    lower = estimates[i],
                    upper = estimates[i + 1L],
                    maximum = maximum,
                    ...
                )
                c(opt[[1L]], opt[[2L]], status)
            },
            double(3L)
        )
    )
    colnames(out) <- c("x", "y", "status")
    out
}

################################################################################
# Helper functions to determine which intervals to search for roots            #
################################################################################

get_search_interval <- function(x_max, y_max, x_min, y_min) {
    # find x-vals where gamma is negative
    neg_min_x <- x_min[y_min < 0]
    # for each of these x-vals, find the closest theta i, where
    # f(theta[i]) > 0 and theta[i] < x-val. Also find theta j,
    # where f(theta[j]) > 0 and theta[j] > x-val.
    interval <- vapply(
        neg_min_x,
        find_closest_thetas,
        x_max = x_max,
        y_max = y_max,
        FUN.VALUE = double(3L)
    )
    # remove duplicate intervals
    keep <- !duplicated(interval[1:2, , drop = FALSE], MARGIN = 2L)
    interval[, keep, drop = FALSE]
}

# This function finds the closest smaller and larger
# value to `minimum` in `x_max` where y_max is positive
find_closest_thetas <- function(minimum, x_max, y_max) {
    cond1 <- y_max > 0
    cond2 <- x_max < minimum
    lower <- max(which(cond1 & cond2))
    cond3 <- x_max > minimum
    upper <- min(which(cond1 & cond3))
    c(
        "lower" = x_max[lower],
        "upper" = x_max[upper],
        "minimum" = minimum
    )
}

################################################################################
# Helper functions that return the lower and upper bound of the                #
# confidence set                                                               #
################################################################################

#' @importFrom stats uniroot
find_lower <- function(estimates_min, SEs_min, f) {
    lower <- estimates_min - SEs_min
    while (f(lower) > 0) {
        lower <- lower - SEs_min
    }
    stats::uniroot(
        f = f,
        lower = lower,
        upper = estimates_min
    )$root
}

#' @importFrom stats uniroot
find_upper <- function(estimates_max, SEs_max, f) {
    upper <- estimates_max + SEs_max
    while (f(upper) > 0) {
        upper <- upper + SEs_max
    }
    stats::uniroot(
        f = f,
        lower = estimates_max,
        upper = upper
    )$root
}

################################################################################
# Helper function that returns a function to optimize                          #
################################################################################

make_function <- function(
    estimates,
    SEs,
    alpha,
    p_fun
) {

    # Add/Overwrite estimate and se args
    function(limit) {
        do.call(
            "p_fun",
            append(
                list(estimates = estimates, SEs = SEs),
                alist(mu = limit)
            )
        ) - alpha
    }
}
