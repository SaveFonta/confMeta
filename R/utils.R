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
    if (!is.na(distr))
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

## Checks thetahat
check_thetahat_arg <- function(thetahat) {
    if (!is_num_fin(thetahat))
        stop(
            "Argument 'thetahat' must be a numeric vector with finite elements."
        )
}

## Checks se
check_se_arg <- function(se, l_thetahat) {
    if (!is_num_fin(se))
        stop("Argument 'se' must be a numeric vector with finite elements.")
    if (min(se) <= 0)
        stop("All entries of argument 'se' must be positive.")
    if (length(se) != l_thetahat && length(se) != 1L)
        stop("Argument 'se' must have length of either 1 or length(thetahat).")
}

## Checks mu
check_mu_arg <- function(mu) {
    if (!is.numeric(mu) || any(!is.finite(mu)) || length(mu) < 1L)
        stop(
            paste0(
                "Argument 'mu' must be a numeric vector of positive length",
                " with finite elements."
            )
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
            )
        )
}

## Checks phi and tau2
check_phiTau2_arg <- function(heterogeneity, phi, tau2) {
    if (heterogeneity == "none") {
        if (!is.null(phi) || !is.null(tau2))
            warning(
                "Ignoring parameter(s) phi and tau2 as heterogeneity = 'none'."
            )
    } else if (heterogeneity == "additive") {
        if (is.null(tau2))
            stop("If heterogeneity = 'additive', tau2 must be provided.")
        if (length(tau2) != 1L || !is_num_fin(tau2))
            stop("Argument 'tau2' must be numeric, finite, and of length 1.")
        if (!is.null(phi))
            warning("Ignoring argument 'phi' as heterogeneity = 'additive'.")
    } else {
        if (is.null(phi))
            stop("If heterogeneity = 'multiplicative', phi must be provided.")
        if (length(phi) != 1L || !is_num_fin(phi) || phi < 0)
            stop("Argument 'phi' must be numeric, finite, and of length 1.")
        if (!is.null(tau2))
            warning(
                "Ignoring argument 'tau2' as heterogeneity = 'multiplicative'."
            )
    }
}

## Checks the alternative argument
check_alternative_arg_hmean <- function(alternative) {
    if (
        length(alternative) != 1L ||
        !(alternative %in% c("none", "less", "greater", "two.sided"))
    )
        stop(
            paste0(
                "Argument 'alternative' must be one of ",
                "c('none', 'less', 'greater', 'two.sided')."
            )
        )
}

## Checks the alternative argument
check_alternative_arg_edg <- function(alternative) {
    if (
        length(alternative) != 1L ||
        !(alternative %in% c("one.sided", "two.sided"))
    )
        stop(
            paste0(
                "Argument 'alternative' must be one of ",
                "c('one.sided', 'two.sided')."
            )
        )
}

## Checks the alternative argument
check_alternative_arg_pearson <- function(alternative) {
    if (
        length(alternative) != 1L ||
        !(alternative %in% c("none"))
    )
        stop(
            paste0(
                "Argument 'alternative' must be one of ",
                "c('none')."
            )
        )
}

## Checks the alternative argument
check_alternative_arg_ktr <- function(alternative) {
    if (
        length(alternative) != 1L ||
        !(alternative %in% c("none"))
    )
        stop(
            paste0(
                "Argument 'alternative' must be one of ",
                "c('none')."
            )
        )
}

## Check the distribution argument used in hMeanChiSqMu()
## - hMeanChiSqMu
check_distr_arg <- function(distr) {
    if (length(distr) != 1L || !(distr %in% c("f", "chisq")))
        stop("Argument 'distr' must be one of c('f', 'chisq').")
}

## Check the w argument used in hMeanChiSqMu()
## - hMeanChiSqMu
check_w_arg <- function(w, thetahat) {
    if (!is_num_fin(w) || length(w) != length(thetahat) || min(w) < 0)
        stop(
            paste0(
                "Argument 'w' must be a numeric vector of the same length as ",
                "'thetahat' with finite and positive elements."
            )
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
    thetahat,
    se,
    heterogeneity,
    phi,
    tau2,
    mu
) {

    # Check thetahat and se are numeric and finite
    ## thetahat
    check_thetahat_arg(thetahat = thetahat)

    ## se
    check_se_arg(se = se, l_thetahat = length(thetahat))

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
        stop("Argument 'level' must be numeric, finite and of length 1.")
    if (level <= 0 || level >= 1)
        stop("Argument 'level' must be between 0 and 1.")
}

# Check wGamma
check_wGamma_arg <- function(wGamma, thetahat) {
    if (length(wGamma) != length(unique(thetahat)) - 1L)
        stop(
            "Argument 'wGamma' must have length length(unique(thetahat)) - 1L."
        )
    if (!is_num_fin(wGamma))
        stop(
            "Argument 'wGamma' must be numeric and all entries must be finite."
        )
}

# Check pValueFUN
check_pValueFUN_arg <- function(pValueFUN) {
    if (!is.function(pValueFUN))
        stop("Argument 'pValueFUN' must be a function.")
}

check_pValueFUNArgs_arg <- function(pValueFUN_args, pValueFUN) {
    if (!is.list(pValueFUN_args))
        stop("Argument 'pValueFUN_args' must be a list'.")
    if ("" %in% names(pValueFUN_args))
        stop("Arument 'pValueFUN_args' must be a named list.")
    # Try to find out what pValueFUN is and check the formals such that the
    # names of the list can be checked
    # However, the above only works for our pValueFUNs since we know their
    # arguments. Custom pValueFUNs might have entirely different Arguments.
}

check_alternative_arg_CI <- function(alternative) {
    if (
        length(alternative) != 1L ||
        !(alternative %in% c("none", "two.sided", "one.sided", "less", "greater"))
    )
        stop(
            paste0(
                "Argument 'alternative' must be one of ",
                "c('none', 'two.sided', 'one.sided', 'less', 'greater')."
            )
        )
}

# Summarise the above functions into one
check_inputs_CI <- function(
    thetahat,
    se,
    level,
    alternative,
    wGamma,
    check_inputs,
    pValueFUN,
    pValueFUN_args
) {
    check_thetahat_arg(thetahat)
    check_se_arg(se = se, l_thetahat = length(thetahat))
    check_level_arg(level = level)
    check_alternative_arg_CI(alternative = alternative)
    check_wGamma_arg(wGamma = wGamma, thetahat = thetahat)
    check_pValueFUN_arg(pValueFUN = pValueFUN)
    check_pValueFUNArgs_arg(pValueFUN_args = pValueFUN_args)
}

################################################################################
# Adjustment of standard errors based on heterogeneity model                   #
################################################################################
# This function adjusts the standard errors depending on the heterogeneity model
# this function is used in the p-value functions:
# - kTRMu
# - hMeanChiSqMu
# - pPearsonMu
adjust_se <- function(se, heterogeneity, phi, tau2) {
    if (heterogeneity == "none")
        se
    else
        switch(
            heterogeneity,
            "additive" = sqrt(se^2 + tau2),
            "multiplicative" = se * sqrt(phi)
        )
}


################################################################################
# Compute the z-values based on thetahat, se and vectorize over mu             #
################################################################################
# This function calculates the z values of thetahat and se for every value of mu
# this function is used in the p-value functions:
# - kTRMu
# - hMeanChiSqMu
# - pPearsonMu
# - pEdgingtonMu

get_z <- function(thetahat, se, mu) {
    n <- length(thetahat)
    z <- vapply(
        mu,
        function(mu) (thetahat - mu) / se,
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

utils::globalVariables(
    c(
        # whereever
        "limit",
        # ggPvalueFunction
        "x", "y", "group", "x_gamma", "y_gamma", "xmin", "xmax", "ymin", "ymax",
        # ForestPlot
        "lower", "upper", "estimate", "id", "color", "name"
    )
)
