get_ci <- function(
    estimates,
    SEs,
    w = NULL,   # [MODIFICA] nuovo argomento per i pesi
    conf_level,
    p_fun
) {
  
  # Keep a copy of estimates and SEs for later (AUCC calculation)
  orig_est <- estimates
  orig_se <- SEs
  orig_w   <- w              # [MODIFICA]<-- aggiungi questa
  
  # f(mu) = p(mu) - alpha
  alpha <- 1 - conf_level
  f <- make_function(
    estimates = estimates,
    SEs = SEs,
    w = w,   # [MODIFICA] passo i pesi
    alpha = alpha,
    p_fun = p_fun
  )
  
  # remove duplicates and sort estimates, SEs, and w
  keep <- !duplicated(estimates)
  estimates <- estimates[keep]
  SEs <- SEs[keep]
  if (!is.null(w)) w <- w[keep]   # [MODIFICA]
  
  o <- order(estimates, decreasing = FALSE)
  estimates <- estimates[o]
  SEs <- SEs[o]
  if (!is.null(w)) w <- w[o]      # [MODIFICA]
  
  # Check if CI exists
  estimates <- matrix(
    c(estimates, f(estimates), rep(0, length(estimates))),
    ncol = 3L,
    dimnames = list(NULL, c("x", "y", "status"))
  )
  
  maxima <- find_optima(estimates = estimates[, 1L], f = f, maximum = TRUE)
  isRelevant_max <- is_relevant(
    f_estimates = estimates[, 2L],
    f_extremum = maxima[, 2L],
    maximum = TRUE
  )
  
  if (any(isRelevant_max)) {
    estimates <- rbind(estimates, maxima[isRelevant_max, ])
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
  
  # Calculate AUCC and ratio
  if (nrow(p_max) > 1L) {
    warning("More than one maximum of the p-value function found. The AUCC ratio is thus undefined.")
    aucc <- NA_real_
    aucc_ratio <- NA_real_
  } else {
    a <- integrate_f(
      max_iter = 7,
      p_fun,
      estimates = orig_est,
      SEs = orig_se,
      w = orig_w,   # [MODIFICA]
      lower = -Inf,
      upper = p_max[, "x"],
      subdivisions = 1000L
    )$value
    b <- integrate_f(
      max_iter = 7,
      p_fun,
      estimates = orig_est,
      SEs = orig_se,
      w = orig_w,   # [MODIFICA]
      lower = p_max[, "x"],
      upper = Inf,
      subdivisions = 1000L
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
    # No CI exists
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
    )
    colnames(out$CI) <- c("lower", "upper")
    colnames(out$gamma) <- c("x", "y")
  } else {
    # CI exists
    estimates_pos <- which(f_estimates > 0)
    idx_min <- min(estimates_pos)
    idx_max <- max(estimates_pos)
    estimates_min <- estimates[idx_min]
    estimates_max <- estimates[idx_max]
    step <- max(SEs)
    
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
    
    estimates <- estimates[idx_min:idx_max]
    f_estimates <- f_estimates[idx_min:idx_max]
    n_intervals <- length(estimates) - 1L
    
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
    
    minima <- gam[, 2L]
    one_pos_theta_only <- length(minima) == 1L && is.na(minima)
    exist_neg_minima <- any(minima < 0)
    search_roots <- !one_pos_theta_only && exist_neg_minima
    if (search_roots) {
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
    
    if (!one_pos_theta_only) {
      gam[, 2L] <- gam[, 2L] + alpha
    }
    
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
