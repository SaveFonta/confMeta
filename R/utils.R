################################################################################
# Wrapper around integrate                                                     #
################################################################################
#' @importFrom stats integrate
#' @noRd

integrate_f <- function (max_iter, ...) {
  exponent <- 0.25
  counter  <- 0L
  
  while (exponent > 0.075 && counter < max_iter) {
    rel_tol <- .Machine$double.eps^exponent
    
    out <- tryCatch(
      integrate(..., rel.tol = rel_tol),
      error = function(e) {
        msg <- conditionMessage(e)
        
        # error due to edgington_w cannot approx properly. To this moment p_edg_w is the only function that can give stop errors
        if (grepl("Exact method infeasible for n > 18", msg, fixed = TRUE)) {
          stop(e)
        }
        
        # failed integration
        structure(list(error = e), class = "retry_integrate")
      }
    )
    
    # not failed integration
    if (!inherits(out, "retry_integrate")) {
      return(out)
    }
    
    # failed integration, retry integration
    exponent <- exponent - 0.025
    counter  <- counter + 1L
  }
  
  stop("integrate_f: integration failed after trying all tolerances.")
}



################################################################################
# Helper functions                                                             #
################################################################################

# Set up a grid of p-value functions with the given
# parameters. This is used in:
# - ggPvalueFunction
# - ForestPlot
make_grid <- function(pValueFUN, heterogeneity, distr) {
    # For each P-value function, make a function
    # that returns a grid of the desired arguments
    # Note: all of these functions must have the same
    # arguments and return a data.frame with the same
    # columns.
    make_grid_hMean <- function(heterogeneity, distr) {
        expand.grid(
            fun_name = "hMeanChiSqMu",
            heterogeneity = heterogeneity,
            distr = distr,
            stringsAsFactors = FALSE
        )
    }
    make_grid_kTRMu <- function(heterogeneity, distr) {
        distr <- NA_character_
        expand.grid(
            fun_name = "kTRMu",
            heterogeneity = heterogeneity,
            distr = distr,
            stringsAsFactors = FALSE
        )
    }
    make_grid_pearson <- function(heterogeneity, distr) {
        distr <- NA_character_
        expand.grid(
            fun_name = "pPearsonMu",
            heterogeneity = heterogeneity,
            distr = distr,
            stringsAsFactors = FALSE
        )
    }
    make_grid_edgington <- function(heterogeneity, distr) {
        distr <- NA_character_
        expand.grid(
            fun_name = "pEdgingtonMu",
            heterogeneity = heterogeneity,
            distr = distr,
            stringsAsFactors = FALSE
        )
    }
    make_grid_fisher <- function(heterogeneity, distr) {
        distr <- NA_character_
        expand.grid(
            fun_name = "pFisherMu",
            heterogeneity = heterogeneity,
            distr = distr,
            stringsAsFactors = FALSE
        )
    }
    # Put functions in a named list
    grid_funs <- list(
        "hMean" = make_grid_hMean,
        "k-Trials" = make_grid_kTRMu,
        "Pearson" = make_grid_pearson,
        "Edgington" = make_grid_edgington,
        "Fisher" = make_grid_fisher
    )
    # Which P-value functions should be included in the grid
    include_funs <- c("hMean", "k-Trials", "Pearson", "Edgington", "Fisher")
    include_funs <- include_funs[include_funs %in% pValueFUN]
    # Subset the list with function to include only those
    # p-value functions that were requested
    grid_funs <- grid_funs[include_funs]
    # Make the grid
    grid <- do.call(
        "rbind",
        lapply(
            seq_along(grid_funs),
            function(i, heterogeneity, distr) {
                fun <- grid_funs[[i]]
                fun_char <- include_funs[i]
                g <- fun(
                    heterogeneity = heterogeneity,
                    distr = distr
                )
                g$pretty_name <- include_funs[i]
                g$name <- make_names(
                    FUN = fun_char,
                    heterogeneity = g$heterogeneity,
                    distr = g$distr
                )
                g
            },
            heterogeneity = heterogeneity,
            distr = distr
        )
    )

    grid
}

# Construct names out of the function and its argument.
# This is used to construct labels for plots etc.
make_names <- function(FUN, heterogeneity, distr) {
    # handle heterogeneity
    heterogeneity <- vapply(
        heterogeneity,
        function(x) {
            switch(
                x,
                "none" = " none",
                "additive" = " add.",
                "multiplicative" = " mult."
            )
        },
        character(1L)
    )
    # handle distr
    do_distr <- !(length(distr) == 1L && is.na(distr))
    if (do_distr)
        distr <- ifelse(is.na(distr), "", paste0(" (", distr, ")"))
    else
        distr <- rep("", length(FUN))
    # Make names
    paste0(FUN, heterogeneity, distr)
}


################################################################################
# Check inputs for p-value functions                                           #
################################################################################
# The following functions are used in:
# - pPearsonMu
# - kTRMu
# - hMeanChiSqMu

## Check whether a vector is numeric and finite and non-NULL
is_num_fin <- function(x) {
    is.numeric(x) && length(x) > 0L && all(is.finite(x))
}

## Checks estimates
check_estimates_arg <- function(estimates) {
    if (!is_num_fin(estimates))
        stop(
            "Argument 'estimates' must be a numeric vector with finite elements.",
            call. = FALSE
        )
}

## Checks SEs
check_SEs_arg <- function(SEs, l_estimates) {
    if (!is_num_fin(SEs))
        stop(
            "Argument 'SEs' must be a numeric vector with finite elements.",
            call. = FALSE
        )
    if (min(SEs) <= 0)
        stop(
            "All entries of argument 'SEs' must be positive.",
            call. = FALSE
        )
    if (length(SEs) != l_estimates && length(SEs) != 1L)
        stop(
            "Argument 'SEs' must have length of either 1 or length(estimates).",
            call. = FALSE
        )
}

## Checks mu
check_mu_arg <- function(mu) {
    if (!is.numeric(mu) || any(!is.finite(mu)) || length(mu) < 1L)
        stop(
            paste0(
                "Argument 'mu' must be a numeric vector of positive length",
                " with finite elements."
            ),
            call. = FALSE
        )
}

## Checks heterogeneity
check_heterogeneity_arg <- function(heterogeneity) {
    if (
        is.null(heterogeneity) ||
        !(heterogeneity %in% c("none", "additive", "multiplicative"))
    )
        stop(
            paste0(
                "Argument 'heterogeneity' must be one of ",
                "c('none', 'additive', 'multiplicative')."
            ),
            call. = FALSE
        )
}

## Checks phi and tau2
check_phiTau2_arg <- function(heterogeneity, phi, tau2) {
    if (heterogeneity == "none") {
        if (!is.null(phi) || !is.null(tau2))
            warning(
                "Ignoring parameter(s) phi and tau2 as heterogeneity = 'none'.",
                call. = FALSE
            )
    } else if (heterogeneity == "additive") {
        if (is.null(tau2))
            stop(
                "If heterogeneity = 'additive', tau2 must be provided.",
                call. = FALSE
            )
        if (length(tau2) != 1L || !is_num_fin(tau2))
            stop(
                "Argument 'tau2' must be numeric, finite, and of length 1.",
                call. = FALSE
            )
        if (!is.null(phi))
            warning(
                "Ignoring argument 'phi' as heterogeneity = 'additive'.",
                call. = FALSE
            )
    } else {
        if (is.null(phi))
            stop(
                "If heterogeneity = 'multiplicative', phi must be provided.",
                call. = FALSE
            )
        if (length(phi) != 1L || !is_num_fin(phi) || phi < 0)
            stop(
                "Argument 'phi' must be numeric, finite, and of length 1.",
                call. = FALSE
            )
        if (!is.null(tau2))
            warning(
                "Ignoring argument 'tau2' as heterogeneity = 'multiplicative'.",
                call. = FALSE
            )
    }
}



# This is is the unified output_p (once called alternative) checker (from Saverio):
check_output_p_arg <- function(output_p) {
  if (
    length(output_p) != 1L ||
    !(output_p %in% c("one.sided", "two.sided"))
  ) {
    stop(
      paste0(
        "Argument 'output_p' must be one of ",
        "c('one.sided', 'two.sided')."
      ),
      call. = FALSE
    )
  }
}

## Check the distribution argument used in hMeanChiSqMu()
## - hMeanChiSqMu
check_distr_arg <- function(distr) {
    if (length(distr) != 1L || !(distr %in% c("f", "chisq")))
        stop("Argument 'distr' must be one of c('f', 'chisq').", call. = FALSE)
}

## Check the w argument used in hMeanChiSqMu()
## - hMeanChiSqMu
check_w_arg <- function(w, estimates) {
    if (!is_num_fin(w) || length(w) != length(estimates) || min(w) < 0)
        stop(
            paste0(
                "Argument 'w' must be a numeric vector of the ",
                "same length as argument ",
                "'estimates' with finite and positive elements."
            ),
            call. = FALSE
        )
}

check_checkInputs_arg <- function(check_inputs) {
    isTRUE(check_inputs) || isFALSE(check_inputs)
}

## Summarise the before functions
## These arguments are present in all p-value functions
# - pPearsonMu
# - hMeanChiSqMu
# - kTRMu
check_inputs_p_value <- function(
    estimates,
    SEs,
    heterogeneity,
    phi,
    tau2,
    mu
) {

    # Check estimates and SEs are numeric and finite
    ## estimates
    check_estimates_arg(estimates = estimates)

    ## SEs
    check_SEs_arg(SEs = SEs, l_estimates = length(estimates))

    # Check mu
    check_mu_arg(mu = mu)

    # Check heterogeneity
    check_heterogeneity_arg(heterogeneity = heterogeneity)

    # Check phi and tau2
    check_phiTau2_arg(heterogeneity = heterogeneity, phi = phi, tau2 = tau2)
}

################################################################################
# Argument checks for hMeanChiSqCI                                             #
################################################################################

# Check level
check_level_arg <- function(level) {
    if (!is_num_fin(level) || length(level) != 1L)
        stop(
            "Argument 'level' must be numeric, finite and of length 1.",
            call. = FALSE
        )
    if (level <= 0 || level >= 1)
        stop(
            "Argument 'level' must be between 0 and 1.",
            call. = FALSE
        )
}

# Check wGamma
# check_wGamma_arg <- function(wGamma, thetahat) {
#     if (length(wGamma) != length(unique(thetahat)) - 1L)
#         stop(
#             "Argument 'wGamma' must have length length(unique(thetahat)) - 1L."
#         )
#     if (!is_num_fin(wGamma))
#         stop(
#             "Argument 'wGamma' must be numeric and all entries must be finite."
#         )
# }

# Check pValueFUN
# check_pValueFUN_arg <- function(pValueFUN) {
#     if (!is.function(pValueFUN))
#         stop("Argument 'pValueFUN' must be a function.")
# }

# check_pValueFUNArgs_arg <- function(pValueFUN_args, pValueFUN) {
#     if (!is.list(pValueFUN_args))
#         stop("Argument 'pValueFUN_args' must be a list'.")
#     if ("" %in% names(pValueFUN_args))
#         stop("Arument 'pValueFUN_args' must be a named list.")
#     # Try to find out what pValueFUN is and check the formals such that the
#     # names of the list can be checked
#     # However, the above only works for our pValueFUNs since we know their
#     # arguments. Custom pValueFUNs might have entirely different Arguments.
# }

check_alternative_arg_CI <- function(alternative) {
    if (
        length(alternative) != 1L ||
        !(alternative %in% c("none", "two.sided", "one.sided", "less", "greater"))
    )
        stop(
            paste0(
                "Argument 'alternative' must be one of ",
                "c('none', 'two.sided', 'one.sided', 'less', 'greater')."
            ),
            call. = FALSE
        )
}

# Summarise the above functions into one
# check_inputs_CI <- function(
#     estimates,
#     SEs,
#     level,
#     alternative,
#     check_inputs,
#     pValueFUN,
#     pValueFUN_args
# ) {
#     check_estimates_arg(estimates)
#     check_SEs_arg(SEs = SEs, l_estimates = length(estimates))
#     check_level_arg(level = level)
#     check_alternative_arg_CI(alternative = alternative)
#     check_pValueFUN_arg(pValueFUN = pValueFUN)
#     check_pValueFUNArgs_arg(pValueFUN_args = pValueFUN_args)
# }

################################################################################
# Adjustment of standard errors based on heterogeneity model                   #
################################################################################
# This function adjusts the standard errors depending on the heterogeneity model
# this function is used in the p-value functions:
# - kTRMu
# - hMeanChiSqMu
# - pPearsonMu
adjust_se <- function(SEs, heterogeneity, phi, tau2) {
    if (heterogeneity == "none")
        SEs
    else
        switch(
            heterogeneity,
            "additive" = sqrt(SEs^2 + tau2),
            "multiplicative" = SEs * sqrt(phi)
        )
}


################################################################################
# Compute the z-values based on estimates, SEs and vectorize over mu           #
################################################################################
# This function calculates the z values of estimates and SEs for every value of
# mu.
# it gives an output a matrix with n rows, where n is number of studies (length(estimates) == length(SEs))
# and it has as many columns as the number of true mu evaluated
# So each row are the Z values of a studies under different true mu.

get_z <- function(estimates, SEs, mu) {
    n <- length(estimates)
    z <- vapply(
        mu,
        function(mu) (estimates - mu) / SEs,
        double(n)
    )
    if (is.null(dim(z))) dim(z) <- c(1L, n)
    z
}

################################################################################
# Global variables                                                             #
# This section is necessary because some of the functions in the ggplot2       #
# package use non-standard evaluation (NSE) which leads to warnings/notes in   #
# R CMD check. Thus we declare all of the variables here.                      #
# ##############################################################################

#' @importFrom utils globalVariables
utils::globalVariables(
    c(
        # autoplot.confMeta
        "x", "y", "group", "xmin", "xmax", "ymin", "ymax",
        "study", "xlim", "lower", "upper", "estimate", "id", "color", "name",
        "y0", "conf_level", "estimates"
    )
)











#=====================================================================================================
# body_p_value_function
#
# This blocks is used in 5 p_value functions:
#        - Tippett
#        - Fisher
#        - Wilkinson
#        - Pearson
#        - Edgington
#       - Edgington w

# the other 2 (Stouffer and Hmean) are based on Z values, 
# so this block cannot be directly used
#=====================================================================================================
#' @noRd
body_p_value_fun <- function(estimates,
                             SEs,
                             mu = 0,
                             heterogeneity = "none",
                             phi = NULL,
                             tau2 = NULL,
                             check_inputs = TRUE,
                             input_p = "greater",
                             output_p = "two.sided",
                             k = rep(1, length(estimates)) 
){
  
  # check inputs
  if (check_inputs) {
    check_inputs_p_value(
      estimates = estimates,
      SEs = SEs,
      mu = mu,
      heterogeneity = heterogeneity,
      phi = phi,
      tau2 = tau2
    )
    check_output_p_arg(output_p = output_p)
    check_k_arg(k = k, estimates = estimates)
  }
  
  # recycle `se` if needed
  if (length(SEs) == 1L) SEs <- rep(SEs, length(estimates))
  
  # adjust se based on heterogeneity model
  SEs <- adjust_se(
    SEs = SEs,
    heterogeneity = heterogeneity,
    phi = phi,
    tau2 = tau2
  )
  
  # get the z-values
  z <- get_z(estimates = estimates, SEs = SEs, mu = mu)
  
  # convert them to p-values
  p <- switch(input_p,
              "two.sided" = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
              "greater"   = stats::pnorm(z, lower.tail = FALSE),
              "less"      = stats::pnorm(z, lower.tail = TRUE),
              stop("input_p must be 'greater','less','two.sided'")
  )
  
  p <- as.matrix(p)
  
  ## adjust for best of k
  if (!all(k == 1L)) {
    p <- adjust_for_best_k(p = p, k = k)
  }
  
  return(p)
}








###################################################################################
##   Best out of k                                                              ###
###################################################################################

# ------------------------------------------------------------------------------
# Input check for k 
# ------------------------------------------------------------------------------

#' @noRd
check_k_arg <- function(k, estimates){
  if (!is.numeric(k) || length(k) != length(estimates) || any(!is.finite(k)) || any(k < 0) )   {
    stop(
      "Argument 'k' must be a numeric vector of the same length of 'estimates', with all elements finite and positive."
    )
  }
  invisible(NULL)
}

# ------------------------------------------------------------------------------
# We are assuming that the observed p value is the smaller out of k studies (with same se)
# The CDF of the minimum of uniform is a Beta: F_k(p_min) = 1 - (1 - p_min)^k 
#
# So if we do p_adj = F_k(p_min) we get that p_adj is Uniform under the null
#
#
# we need to adjust rowise. since each column is a mu value and
# row 1 = study 1 needs to be adjusted with k[1], row 2 with k[2] etc...
#
# R works perfectly. Try to believe: 
#
## m <- matrix(rep(2,9), nrow = 3, ncol = 3)
## k <- c(1,2,3)
## (m+1)^k
# ------------------------------------------------------------------------------

#' @noRd
adjust_for_best_k <- function(p, k) {
  #1 - (1 - p)^k
  # log1p(x) computes log(1 + x), while expm1 computes exp(x) - 1, so we do -expm to get 1 - exp(x)
  - expm1(k * log1p(-p))
}







# ------------------------------------------------------------------------------
# best-out-of-k for methods that use Z scores and not p values
#
# The idea is z->p->adjust->z.
# The formula depends on input_p because the direction determines which
# tail of the normal distribution the one-sided p-value uses
#
# Derivations:
#
#   "greater": p_i = 1 - Phi(z_i)
#              p_adj = 1 - (1 - p_i)^k = 1 - (1 - 1 - Phi(z_i))^k = 1 - Phi(z_i)^k
#              z_adj = Phi^{-1}(1 - p_adj) = (replace p_adj) =  Phi^{-1}(Phi(z_i)^k)
#
#   "less":    p_i = Phi(z_i)
#              p_adj = 1 - (1 - p_i)^k = 1 - (1 - Phi(z_i))^k = 1 - (Phi(-z_i))^k 
#              z_adj = Phi^{-1}(p_adj) = (replace p_adj) = Phi^{-1}(1 - (Phi(-z_i))^k =  -Phi^{-1}(Phi(-z_i)^k)
#
#   "two.sided": p_i = 2*Phi(-|z_i|)
#              p_adj = 1 - (1 - 2*Phi(-|z_i|))^k
#              z_adj = Phi^{-1}(1 - p_adj) * sign(z_i)
# not sure baout two sided, but not relevant actually
# ------------------------------------------------------------------------------

#' @importFrom stats pnorm qnorm
#' @noRd
adjust_z_for_best_k <- function(z, k, input_p) {
  z <- as.matrix(z)
  switch(
    input_p,
    "greater" = {
      #stats::qnorm(stats::pnorm(z)^k)
      stats::qnorm(k * stats::pnorm(z, log.p = TRUE), log.p = TRUE)
    },
    "less" = {
      #-stats::qnorm(stats::pnorm(-z)^k)
      -stats::qnorm(k * stats::pnorm(-z, log.p = TRUE), log.p = TRUE)
    },
    "two.sided" = { #didn't use log space, but no function that combines z values can be two sided for now, so we never get to this block 
            p_two <- 2 * stats::pnorm(-abs(z))
            p_adj <- 1 - (1 - p_two)^k
            stats::qnorm(1 - p_adj) * sign(z)
    },
    stop("input_p must be 'greater', 'less', or 'two.sided'")
  )
}








# ------------------------------------------------------------------------------
# Returns a single character string: "greater", "less", or "two.sided".
# ------------------------------------------------------------------------------

#' @noRd
resolve_input_p <- function(fun, ell) {
  
  
  # note that since we use remove_unused, if the p_function doesnt have a input_p parameter,
  # even if the user specify it in the confMeta call, it gets deleted, and it will fall to the other cases 

  # if input_p is specified by user
  
  if ("input_p" %in% names(ell)) {
    return(ell[["input_p"]])
  }
  
  # fun's own default for input_p, if it has one
  f <- formals(fun)
  if ("input_p" %in% names(f)) {
    default_val <- f[["input_p"]]

    resolved <- tryCatch(eval(default_val), error = function(e) NULL)
    if (!is.null(resolved) && is.character(resolved)) {
      return(resolved[1L])
    }
  }
  
  # fun has no input_p argument (e.g. p_hmean) -> "greater"
  "greater"
}





# ------------------------------------------------------------------------------
# Adjusted individual study summaries
#
#
# for input_p == ("greater" or "less"):
#
#     theta_k(q) = estimates_i - SEs_i * qnorm( (1 - q)^(1/k_i) )
#
#   Evaluated at q = alpha/2, 0.5, 1-alpha/2 gives upper, adjusted_estimate,
#   lower respectively.
#
# NON-DIRECTIONAL selection input_p ==("two.sided"):
# means that the reported experiment is the most extreme result in either direction
# 
#    # The adjusted two-sided p-value function
#   F(theta0) = 1 - (2*Phi(|thetahat - theta0|/se) - 1)^k
# is symmetric and unimodal, so it is not a cdf  
#
# "decentralize" it to the monotone CDF:
#
#   F_tilde(theta0) = a + (-1)^a / 2 * F(theta0)
#
# where a = I(theta0 > thetahat). 
#
# For theta0 < thetahat (lower boundary, a=0):
#   F_tilde = F/2 = alpha/2  =>  F = alpha
#   => (2*Phi(|z|) - 1)^k = 1 - alpha
#   => |z| = qnorm((1 + (1-alpha)^(1/k)) / 2)
#
# By symmetry the upper boundary uses the same z_alpha.
# The point estimate is the median of F_tilde, which is thetahat
# (no shrinkage). For k=1: z_alpha = qnorm(1 - alpha/2). 
# ------------------------------------------------------------------------------


#' @importFrom stats qnorm
#' @noRd
adjusted_individual_cis <- function(estimates, SEs, conf_level, k, input_p) {
  
  alpha <- 1 - conf_level
  
  if (input_p == "two.sided") {
    
    
    
    z_alpha           <- stats::qnorm((1 + (1 - alpha)^(1 / k)) / 2)
    lower             <- estimates - SEs * z_alpha
    adjusted_estimate <- estimates
    upper             <- estimates + SEs * z_alpha
    
  } else if (input_p == "greater") {
    
    #   theta_k(q) = estimates - SEs * qnorm( (1-q)^(1/k) )
    #    # point estimate: theta_k(0.5). Shrinks for k > 1.
    lower             <- estimates - SEs * stats::qnorm((1 - alpha / 2)^(1 / k))
    adjusted_estimate <- estimates - SEs * stats::qnorm(0.5^(1 / k))
    upper             <- estimates - SEs * stats::qnorm((alpha / 2)^(1 / k))
    
  } else if (input_p == "less") {
    
    lower             <- estimates + SEs * stats::qnorm((alpha / 2)^(1 / k))
    adjusted_estimate <- estimates + SEs * stats::qnorm(0.5^(1 / k))
    upper             <- estimates + SEs * stats::qnorm((1 - alpha / 2)^(1 / k))
    
  } else {
    stop("input_p must be 'greater', 'less', or 'two.sided'", call. = FALSE)
  }
  
  matrix(
    c(lower, adjusted_estimate, upper),
    ncol = 3L,
    dimnames = list(
      names(estimates),
      c("lower", "adjusted_estimate", "upper")
    )
  )
}  















# ------------------------------------------------------------------------------
# Since the adjusted estimates in best of k are not Normal (skewed CIs). We need to compute the mean and the 
# se of the individual studies (adjusted for best of k), to use them to compute reference methods (FE/RE/HK/HC) 

# The idea is to obtain the confidence density of each study (derivative of one sided P val function) and find 
# mean and se of each of tose
#
# The confidence density for each study is:
#
#   "greater":
#     f(theta0) = k * Phi((thetahat - theta0)/se)^(k-1) * phi((thetahat - theta0)/se) / se
#
#   "less" (mirror of "greater"):
#     f(theta0) = k * Phi((theta0 - thetahat)/se)^(k-1) * phi((theta0 - thetahat)/se) / se
#
#   "two.sided" (density of F_tilde, the decentralized monotone CDF):
#     f(theta0) = k * (2*Phi(|thetahat - theta0|/se) - 1)^(k-1) * phi(|thetahat - theta0|/se) / se

#
# For k = 1 all three reduce exactly to N(thetahat, se^2),
# recovering the standard normal confidence density.
#
# ------------------------------------------------------------------------------

#' @importFrom stats integrate dnorm pnorm
#' @noRd
confidence_density_mean_se <- function(estimates, SEs, k, input_p) {
  
  n <- length(estimates)
  
  out <- matrix(
    NA_real_,
    nrow = n,
    ncol = 2L,
    dimnames = list(names(estimates), c("mean", "se"))
  )
  
  for (i in seq_len(n)) {
    
    thetahat <- estimates[i]
    se       <- SEs[i]
    ki       <- k[i]
    
    # Integration bounds: +-20 SEs
    lo <- thetahat - 20 * se
    hi <- thetahat + 20 * se
    
    # Confidence density for study i
    f <- switch(
      input_p,
      "greater" = function(theta0) {
        z <- (thetahat - theta0) / se
        ki * stats::pnorm(z)^(ki - 1) * stats::dnorm(z) / se
      },
      "less" = function(theta0) {
        z <- (theta0 - thetahat) / se
        ki * stats::pnorm(z)^(ki - 1) * stats::dnorm(z) / se
      },
      "two.sided" = function(theta0) {
        z <- abs(thetahat - theta0) / se
        ki * (2 * stats::pnorm(z) - 1)^(ki - 1) * stats::dnorm(z) / se
      },
      stop("input_p must be 'greater', 'less', or 'two.sided'", call. = FALSE)
    )
    
    # Mean: integral of theta0 * f(theta0)
    mu_i <- integrate_f(
      max_iter = 10L,
      f    = function(theta0) theta0 * f(theta0),
      lower = lo,
      upper = hi
    )$value
    
    # Variance: integral of (theta0 - mu)^2 * f(theta0)
    var_i <- integrate_f(
      max_iter = 10L,
      f    = function(theta0) (theta0 - mu_i)^2 * f(theta0),
      lower = lo,
      upper = hi
    )$value
    
    out[i, "mean"] <- mu_i
    out[i, "se"]   <- sqrt(var_i)
  }
  
  out
}
