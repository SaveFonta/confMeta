#' @title Creating *confMeta* objects
#' @description Function to create objects of class `confMeta`. This is the
#'     main class within the package. For an overview of available methods
#'     run `methods(confMeta)`.
#' @param estimates A vector containing the normalized individual effect
#'     estimates. Must be of the same length as `SEs` and coercible to
#'     type 'double'.
#' @param SEs The standard errors of the normalized individual effect estimates.
#'     Must be of the same length as `estimates` and coercible to
#'     type 'double'.
#' @param study_names Either `NULL` (the default) or a vector that can
#'     be coerced to type 'character'. Must be of the same length as
#'     arguments `estimates` and `SEs`. The vector is used to
#'     differentiate the individual studies. For more information, check out the
#'     Details section.
#' @param conf_level The confidence level. Must be a numeric vector of length
#'     one with its value being in (0, 1).
#' @param fun A function that combines individual effect estimates and the
#'     corresponding standard errors into a combined p-value. The function must
#'     have arguments named 'estimates', 'SEs', and 'mu'.
#'     Additional arguments are also allowed, but these must either have a
#'     default value or be passed via the `...` argument.
#'     For more information, see the Details section.
#' @param fun_name A character vector of length 1. The vector serves as an
#'     identifier for the function `fun` and is only used as a label in
#'     plots. The default is just the literal that the function `fun` is
#'     bound to.
#' @param w An optional numeric vector of weights for the studies. Must be of 
#'     the same length as `estimates` and `SEs`, finite, and non-negative. 
#'     If `NULL` (the default), equal weights are assumed. Used by weighted 
#'     p-value combination functions such as the weighted Edgington method.
#' @param ... Additional arguments passed to `fun`. See the Details
#'     section.
#' @return An S3 object of class `confMeta`. The object contains the
#'     following elements:
#'     \describe{
#'         \item{estimates}{The normalized individual effect estimates.}
#'         \item{SEs}{The standard errors of the normalized individual effect
#'             estimates.}
#'         \item{study_names}{The names of the individual studies.}
#'         \item{conf_level}{The confidence level.}
#'         \item{w}{The weights for the studies (if specified). If `NULL`, 
#'             equal weights were assumed.}
#'         \item{individual_cis}{The confidence intervals for the
#'             individual effects. The exact calculation of these intervals
#'             can be found in the Details section.}
#'         \item{p_fun}{A function with arguments named 'estimates', 'SEs',
#'             'conf_level', and 'mu'. This is the p-value function that is
#'             used to find the combined confidence intervals.}
#'         \item{fun_name}{The name of the function. It is only used in plots
#'             as a legend entry.}
#'         \item{joint_cis}{The combined confidence interval(s). These are
#'             calculated by finding the mean values where the $p$-value
#'             function is larger than the confidence level in element
#'             `conf_level`.}
#'         \item{gamma}{The local minima within the range of the individual
#'             effect estimates. Column 'x' refers to the mean `mu` and
#'             column 'y' contains the corresponding $p$-value.}
#'         \item{p_max}{The local maxima of the $p$-value function. The
#'             column 'x' refers to the mean `mu` and the column 'y' contains
#'             the corresponding $p$-value.}
#'         \item{p_0}{The value of the $p$-value at `mu` = 0}
#'         \item{comparison_cis}{Combined confidence intervals calculated
#'             with other methods. These can be used for comparison
#'             purposes. Currently, these other methods are random effects
#'             (REML), Hartung & Knapp, and Henmi & Copas.}
#'         \item{comparison_p_0}{The same as in element 'p_0' but for the
#'             comparison methods (Random effects, Hartung & Knapp,
#'             Henmi & Copas).}
#'     }
#' @details
#'     # Function arguments
#'
#'     The argument `study_names` is used to differentiate between the
#'     different individual estimates. If the argument is set to `NULL`,
#'     the element 'study_names' in the return object will just be a
#'     character vector with elements "Study n" where n is a number from 1
#'     to `length(estimates)`. These names are only used in some of the
#'     `autoplot` methods.
#'
#'     The argument `fun` must have arguments 'estimates', 'SEs',
#'     and 'mu' but it can also have further arguments.
#'     However, these must either have a default value or
#'     need to be passed via the `...` argument. If there
#'     are additional arguments passed via `...`,
#'     `confMeta` will internally create a new
#'     function that calls `fun` with the additional arguments fixed.
#'     Thus, any argument passed via `...` overwrites existing
#'     defaults.
#'     Since this is the p-value function that is used to calculate the
#'     combined confidence interval(s), it should return a vector of class
#'     'numeric' with value(s) in the interval [0, 1].
#'
#'     The argument `w` allows assigning weights to the studies. If not 
#'     specified, all studies receive equal weight (`w = rep(1, n)`). 
#'     Weighted p-value combination functions (e.g. weighted Edgington) 
#'     will use these values directly; other functions that do not depend 
#'     on weights will ignore `w`.
#'
#'     # Confidence intervals
#'
#'     The `individual_cis` are calculated as
#'     \deqn{x_{i} \pm \Phi^{-1}{\text{conf_level}} \cdot \sigma_{i}}
#'     where \eqn{x_{i}} corresponds to the elements of vector
#'     `estimates`, \eqn{\Phi^{-1}} is the quantile function of
#'     the standard normal distribution, conf_level is the confidence
#'     level passed as argument `conf_level`, and
#'     \eqn{\sigma_{i}}, are the standard errors passed in argument
#'     `SEs`.
#'
#'     The boundaries of the confidence intervals returned in element
#'     `joint_cis` are found by searching where the function returned
#'     in element 'p_fun' is equal to 1-`conf_level`.
#' 
#' @examples
#'     # Simulate effect estimates and standard errors
#'     set.seed(42)
#'     n <- 5
#'     estimates <- rnorm(n)
#'     SEs <- rgamma(n, 5, 5)
#'     conf_level <- 0.95
#'
#'     # Construct a simple confMeta object using p_edgington as
#'     # the p-value function
#'     cm <- confMeta(
#'         estimates = estimates,
#'         SEs = SEs,
#'         conf_level = conf_level,
#'         fun = p_edgington,
#'         fun_name = "Edgington  (one-sided input)",
#'         input_p = "one.sided"
#'     )
#'     cm2 <- confMeta(
#'         estimates = estimates,
#'         SEs = SEs,
#'         conf_level = conf_level,
#'         fun = p_edgington,
#'         fun_name = "Edgington (two-sided input)",
#'         input_p = "two.sided"
#'     )
#'
#'     # Plot the object
#'     autoplot(cm, cm2, type = "p")                   # p-value function plot
#'     autoplot(cm, cm2, type = "forest")              # forest plot
#'     autoplot(cm, cm2, type = c("p", "forest"))      # both
#'
#' @export
confMeta <- function(
    estimates,
    SEs,
    study_names = NULL,
    conf_level = 0.95,
    fun,
    fun_name = NULL,
    w = NULL,   # [MOD] opzionale, può rimanere NULL
    ...
) {
  
  if (is.null(study_names)) {
    study_names <- paste0("Study ", seq_along(estimates))
  }
  if (is.null(fun_name)) {
    fun_name <- deparse1(substitute(fun))
  }
  
  ell <- list(...)
  ell <- remove_unused(fun = fun, ell = ell)
  
  if (inherits(estimates, "numeric")) estimates <- as.double(estimates)
  if (inherits(SEs, "numeric")) SEs <- as.double(SEs)
  if (inherits(conf_level, "numeric")) conf_level <- as.double(conf_level)
  study_names <- as.character(study_names)
  
  validate_inputs(
    estimates = estimates,
    SEs = SEs,
    study_names = study_names,
    conf_level = conf_level,
    fun = fun,
    fun_name = fun_name,
    ell = ell,
    w = w #[MOD]
  )
  
  p_fun <- make_p_fun(
    fun = fun,
    ell = ell
  )
  
  new_confMeta(
    estimates = estimates,
    SEs = SEs,
    w = w,  # [MOD] può anche essere NULL, ci penserà new_confMeta
    study_names = study_names,
    conf_level = conf_level,
    p_fun = p_fun,
    fun_name = fun_name
  )
}

# Constructor function
#' @importFrom stats qnorm
new_confMeta <- function(
    estimates = double(),
    SEs = double(),
    w = NULL,   # [MOD] 
    study_names = character(),
    conf_level = double(1L),
    p_fun,
    fun_name
) {
  
  
  # Calculate individual CIs (classic Wald type)
  alpha <- 1 - conf_level
  se_term <- stats::qnorm(1 - alpha / 2) * SEs 
  individual_cis <- matrix(
    c(
      estimates - se_term,
      estimates + se_term
    ), 
    ncol = 2L,
    dimnames = list(study_names, c("lower", "upper"))
  )
  
  # Calculate the joint CIs with the p-value function
  joint_cis <- get_ci(
    estimates = estimates,
    SEs = SEs,
    w = w,  # [MOD] se non è NULL, verrà usato dentro get_ci
    conf_level = conf_level,
    p_fun = p_fun
  )
  
  # Calculate joint CIs with the comparison methods
  method <- c("fe", "re", "hk", "hc")
  comparison <- get_stats_others(
    method = method,
    estimates = estimates,
    SEs = SEs,
    conf_level = conf_level
  )
  
  # Return object
  structure(
    list(
      estimates = estimates,
      SEs = SEs,
      w = w,   # [MOD] salvo i pesi se presenti
      study_names = study_names,
      conf_level = conf_level,
      p_fun = p_fun,
      fun_name = fun_name,
      individual_cis = individual_cis,
      joint_cis = joint_cis$CI,
      gamma = joint_cis$gamma,
      p_max = joint_cis$p_max,
      p_0 = joint_cis$p_0,
      aucc = joint_cis$aucc,
      aucc_ratio = joint_cis$aucc_ratio,
      comparison_cis = comparison$CI,
      comparison_p_0 = comparison$p_0
    ),
    class = "confMeta"
  )
}

# Input checker
validate_inputs <- function(
    estimates,
    SEs,
    study_names,
    conf_level,
    fun,
    fun_name,
    ell,
    w = NULL
) {
  if (is.null(w)) {
    check_equal_length(             # estimates, SEs, study_names should
      estimates = estimates,        # have same length
      SEs = SEs,
      study_names = study_names
    )
  } else {         
    check_equal_length(
      estimates = estimates,
      SEs = SEs,
      study_names = study_names,
      w = w
    )
  }
  check_length_1(x = conf_level)         # conf_level must be of length 1
  check_length_1(x = fun_name)           # fun_name must be of length 1
  
  # Check validity of values
  check_all_finite(x = estimates)        # no NAs, NaNs etc in estimates
  check_all_finite(x = SEs)              # no NAs, NaNs etc in SEs
  check_all_finite(x = conf_level)       # no NAs, NaNs etc in conf_level
  if (!is.null(w)) {
    check_all_finite(x = w)             
    if (any(w < 0)) {
      stop("Weights (w) must be non-negative.") #[MOD] 
    }
  }  
  check_prob(x = conf_level)             # conf_level must be between 0 & 1
  check_fun_args(fun = fun, ell = ell)   # function must have correct args
  
  # Check the function and its arguments
  invisible(NULL)
}

# Validator function
validate_confMeta <- function(confMeta) {
  
  # All valid names (aggiunto w)
  cm_elements <- c(
    "estimates",
    "SEs",
    "study_names",
    "conf_level",
    "p_fun",
    "fun_name",
    "individual_cis",
    "joint_cis",
    "gamma",
    "p_max",
    "p_0",
    "aucc",
    "aucc_ratio",
    "comparison_cis",
    "comparison_p_0",
    "w"   # [MOD]
  )
  
  ok <- cm_elements %in% names(confMeta)
  if (!all(ok)) {
    miss <- cm_elements[!ok]
    msg <- paste0(
      "The confMeta object is missing components: ",
      format_elements(miss)
    )
    stop(msg)
  }
  
  # type checks
  with(
    confMeta,
    {
      # Check types
      check_type(x = estimates, "double", val = TRUE)
      check_type(x = SEs, "double", val = TRUE)
      check_type(x = study_names, "character", val = TRUE)
      check_type(x = fun_name, "character", val = TRUE)
      check_type(x = conf_level, "double", val = TRUE)
      check_type(x = individual_cis, "double", val = TRUE)
      check_type(x = joint_cis, "double", val = TRUE)
      check_type(x = gamma, "double", val = TRUE)
      check_type(x = p_max, "double", val = TRUE)
      check_type(x = p_0, "double", val = TRUE)
      check_type(x = aucc, "double", val = TRUE)
      check_type(x = aucc_ratio, "double", val = TRUE)
      check_type(x = comparison_cis, "double", val = TRUE)
      check_type(x = comparison_p_0, "double", val = TRUE)
      if (!is.null(w)) check_type(x = w, "double", val = TRUE)   # [MOD]
      check_is_function(x = p_fun)
      
      # Check classes
      check_class(x = individual_cis, class = "matrix", val = TRUE)
      check_class(x = joint_cis, class = "matrix", val = TRUE)
      check_class(x = gamma, class = "matrix", val = TRUE)
      check_class(x = p_max, class = "matrix", val = TRUE)
      check_class(x = p_0, class = "matrix", val = TRUE)
      check_class(x = aucc, class = "numeric", val = TRUE)
      check_class(x = aucc_ratio, class = "numeric", val = TRUE)
      check_class(x = comparison_cis, "matrix", val = TRUE)
      check_class(x = comparison_p_0, "matrix", val = TRUE)
      
      # Check validity of values
      check_all_finite(x = estimates, val = TRUE)
      check_all_finite(x = SEs, val = TRUE)
      check_all_finite(x = conf_level, val = TRUE)
      if (!is.null(w)) check_all_finite(x = w, val = TRUE)         # [MOD]
      check_prob(x = conf_level, val = TRUE)
      check_fun_args(
        fun = p_fun,
        ell = list(),
        val = TRUE
      )
      
      # Check lengths
      if (is.null(w)) { #[MOD]
        check_equal_length(
          estimates = estimates,
          SEs = SEs,
          study_names = study_names
        )
      } else {
        check_equal_length(
          estimates = estimates,
          SEs = SEs,
          study_names = study_names,
          w = w
        )
      }
      check_length_1(x = conf_level)
      check_length_1(x = fun_name)
      check_length_1(x = aucc)
      check_length_1(x = aucc_ratio)
      
      invisible(NULL)
    }
  )
  
  invisible(NULL)
}

# ==============================================================================
# Calculate the CIs using the other methods ====================================
# ==============================================================================

#' @importFrom meta metagen
get_obj_re <- function(estimates, SEs, conf_level) {
    meta::metagen(
        TE = estimates, seTE = SEs, sm = "MD",
        level = conf_level, method.tau = "REML",
        random = TRUE, common = FALSE
    )
}

#' @importFrom meta metagen
get_obj_fe <- function(estimates, SEs, conf_level) {
    meta::metagen(
        TE = estimates, seTE = SEs, sm = "MD",
        level = conf_level, method.tau = "REML",
        random = FALSE, common = TRUE
    )
}

#' @importFrom meta metagen
get_obj_hk <- function(estimates, SEs, conf_level) {
    meta::metagen(
        TE = estimates, seTE = SEs, sm = "MD",
        level = conf_level, method.tau = "REML", method.random.ci = "HK",  #[MOD]--> hakn = TRUE is deprecated
        common = FALSE, random = TRUE
    )
}

#' @importFrom metafor hc rma
get_obj_hc <- function(estimates, SEs, conf_level) {
    metafor::hc(
        object = metafor::rma(yi = estimates, sei = SEs, level = conf_level)
    )
}

# Get the CI
get_ci_re <- function(re) {
    matrix(
        c(re$lower.random, re$upper.random),
        ncol = 2L,
        dimnames = list(get_method_names()["re"], c("lower", "upper"))
    )
}

get_ci_fe <- function(fe) {
    matrix(
        c(fe$lower.common, fe$upper.common),
        ncol = 2L,
        dimnames = list(get_method_names()["fe"], c("lower", "upper"))
    )
}

get_ci_hk <- function(hk) {
    matrix(
        c(hk$lower.random, hk$upper.random),
        ncol = 2L,
        dimnames = list(get_method_names()["hk"], c("lower", "upper"))
    )
}

get_ci_hc <- function(hc) {
    matrix(
        c(hc$ci.lb, hc$ci.ub),
        ncol = 2L,
        dimnames = list(get_method_names()["hc"], c("lower", "upper"))
    )
}

# Get the p-value for the null-effect
get_pval_re <- function(obj) {
    matrix(
        c(0, obj$pval.random),
        ncol = 2L,
        dimnames = list(get_method_names()["re"], c("x", "y"))
    )
}

get_pval_fe <- function(obj) {
    matrix(
        c(0, obj$pval.common),
        ncol = 2L,
        dimnames = list(get_method_names()["fe"], c("x", "y"))
    )
}

get_pval_hk <- function(obj) {
    matrix(
        c(0, obj$pval.random),
        ncol = 2L,
        dimnames = list(get_method_names()["hk"], c("x", "y"))
    )
}

#' @importFrom ReplicationSuccess ci2p
get_pval_hc <- function(obj, conf_level) {
    ci <- get_ci_hc(obj)
    p <- ReplicationSuccess::ci2p(
        lower = ci[, "lower"],
        upper = ci[, "upper"],
        conf.level = conf_level,
        ratio = FALSE,
        alternative = "two.sided"
    )
    matrix(
        c(0, p),
        ncol = 2L,
        dimnames = list(get_method_names()["hc"], c("x", "y"))
    )
}

get_method_names <- function() {
    c(
        "fe" = "Fixed effect",
        "re" = "Random effects",
        "hk" = "Hartung & Knapp",
        "hc" = "Henmi & Copas"
    )
}

# Calculate everything
get_stats_others <- function(method, estimates, SEs, conf_level) {
    nms <- get_method_names()
    stopifnot(all(method %in% names(nms)))
    names(method) <- method
    res <- lapply(
        method,
        function(x, estimates, SEs, conf_level, nms) {
            obj <- switch(
                x,
                "re" = get_obj_re(
                    estimates = estimates,
                    SEs = SEs,
                    conf_level = conf_level
                ),
                "hk" = get_obj_hk(
                    estimates = estimates,
                    SEs = SEs,
                    conf_level = conf_level
                ),
                "hc" = get_obj_hc(
                    estimates = estimates,
                    SEs = SEs,
                    conf_level = conf_level
                ),
                "fe" = get_obj_fe(
                    estimates = estimates,
                    SEs = SEs,
                    conf_level = conf_level
                )
            )
            ci <- switch(
                x,
                "re" = get_ci_re(re = obj),
                "hk" = get_ci_hk(hk = obj),
                "hc" = get_ci_hc(hc = obj),
                "fe" = get_ci_fe(fe = obj)
            )
            pval <- switch(
                x,
                "re" = get_pval_re(obj = obj),
                "hk" = get_pval_hk(obj = obj),
                "hc" = get_pval_hc(obj = obj, conf_level = conf_level),
                "fe" = get_pval_fe(obj = obj)
            )
            list(
                ci = ci,
                p_0 = pval
            )
        },
        estimates = estimates,
        SEs = SEs,
        conf_level = conf_level,
        nms = nms
    )

    cis <- do.call("rbind", lapply(res, "[[", i = "ci"))
    p <- do.call("rbind", lapply(res, "[[", i = "p_0"))

    list("CI" = cis, "p_0" = p)
}

# ==============================================================================
# Remove unnecessary arguments from ... ========================================
# ==============================================================================

remove_unused <- function(fun, ell) {
    fa <- names(formals(fun))
    ell_args <- names(ell)
    # Check whether there are unused arguments and if so, give a message
    unused_args <- !(ell_args %in% fa)
    if (any(unused_args)) {
        remove <- ell_args[unused_args]
        message(
            paste0(
                "Ignoring unused argument(s) ",
                format_elements(remove),
                "."
            )
        )
        ell <- ell[!unused_args]
    }
    ell
}

# ==============================================================================
# Make the p-value function  ===================================================
# ==============================================================================

make_p_fun <- function(fun, ell) {

    f <- formals(fun)                           # get the arguments of fun

    if (length(ell) > 0L) {                     # in case of non-empty dotargs
        f <- f[!(names(f) %in% names(ell))]     # remove args passed via ...
        f <- append(f, ell)                     # add these args from ...
    }

    out <- fun
    formals(out) <- f

    out
}


# ==============================================================================
# Argument checks ==============================================================
# ==============================================================================

################################################################################
# Checking whether x is a function                                             #
################################################################################

check_is_function <- function(x) {
    if (!is.function(x)) {
        obj <- deparse1(substitute(x))
        msg <- paste0("Argument `", obj, "` must be a function.")
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

################################################################################
# Checking the type of a variable                                              #
################################################################################

check_all_finite <- function(x, val = FALSE) {
    if (!all(is.finite(x))) {
        obj <- deparse1(substitute(x))
        msg <- paste0(
            if (length(x) > 1L) {
                paste0(
                    "All entries of ",
                    if (val) "element `" else "argument `"
                )
            } else {
                if (val) "Element `" else "Argument `"
            },
            obj,
            "` must be finite."
        )
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

################################################################################
# Checking the type of a variable                                              #
################################################################################

check_type <- function(x, type, val = FALSE) {
    if (!(typeof(x) == type)) {
        obj <- deparse1(substitute(x))
        msg <- paste0(
            if (val) "Element `" else "Argument `",
            obj,
            "` must be of type '",
            type,
            "'."
        )
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

################################################################################
# Checking the class of a variable                                             #
################################################################################

check_class <- function(x, class, val = FALSE) {
    if (!inherits(x, class)) {
        obj <- deparse1(substitute(x))
        msg <- paste0(
            if (val) "Element `" else "Argument `",
            obj,
            "` must be of class '",
            class,
            "'."
        )
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

################################################################################
# Check that the function has correct arguments                                #
################################################################################

check_fun_args <- function(fun, ell = NULL, val = FALSE) {

    # fun must have the following arguments
    must_have_args <- c(
        "estimates",
        "SEs",
        "mu"
    )

    # Get info about the function and its arguments
    f <- formals(fun)              # formals(fun)
    fa <- names(f)                 # formalArgs(fun)

    # Get the additional arguments
    if (!is.null(ell)) {
        ell_args <- names(ell)     # names of the ellipsis args
    }

    # Check whether the function has all the required
    ok_required <- has_all(
        x = must_have_args,
        superset = fa
    )
    if (!ok_required) {
        msg <- paste0(
            "Function in ",
            if (val) "element `p_fun` " else "argument `fun` ",
            "must have arguments ",
            "named 'estimates', 'SEs', and 'mu'."
        )
        stop(msg, call. = FALSE)
    }

    # Check the non-required arguments whether they have a default
    # or are passed via ...
    add_args <- fa[!(fa %in% must_have_args)]  # Get additional args
    not_ok <- vapply(                          # Check whether they have default
        f[add_args],
        is_missing,
        logical(1L),
        USE.NAMES = TRUE
    )
    if (any(not_ok)) {                         # If any w/o default, check ...
        if (!val) {
            args_check <- names(not_ok)[not_ok]
            args_ok <- args_check %in% ell_args
            names(args_ok) <- args_check
            not_ok <- !args_ok
        }
        if (any(not_ok)) {
            msg <- paste0(
                "Invalid ",
                if (val) "element `p_fun`" else "argument `fun`",
                ". The function argument(s) ",
                format_elements(names(not_ok)[not_ok]),
                " must ",
                if (val) {
                    "have a default."
                } else {
                    "either have a default or be passed via `...`."
                }
            )
            stop(msg, call. = FALSE)
        }
    }

    invisible(NULL)
}

is_missing <- function(x) {
    if (!all(nzchar(x)) && is.name(x)) {
        TRUE
    } else {
        FALSE
    }
}

has_all <- function(x, superset) {
    all(x %in% superset)
}

################################################################################
# Checking for equal lengths of variables                                      #
################################################################################

check_equal_length <- function(..., val = FALSE) {
    arguments <- list(...)
    l <- vapply(arguments, length, integer(1L), USE.NAMES = TRUE)
    ind <- l == l[1L]
    ok <- all(ind)
    if (!ok) {
        nms <- names(l)
        msg <- paste0(
            if (val) "Elements " else "Arguments ",
            format_elements(nms),
            " must have the same length."
        )
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

format_elements <- function(x) {
    nms <- paste0("`", x, "`")
    last <- length(x)
    if (last == 1L) {
        nms
    } else if (last == 2L) {
        paste0(nms[1L], " and ", nms[2L])
    } else {
        nms_first <- paste0(nms[-last], collapse = ", ")
        nms_last <- paste0(", and ", nms[last])
        paste0(nms_first, nms_last)
    }
}

################################################################################
# Check for length equals 1                                                    #
################################################################################

check_length_1 <- function(x, val = FALSE) {
    if (length(x) != 1L) {
        obj <- deparse1(substitute(x))
        msg <- paste0(
            if (val) "Element `" else "Argument `",
            obj,
            "` must be of length 1."
        )
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

################################################################################
# Check for value between 0 and 1                                              #
################################################################################

check_prob <- function(x, val = FALSE) {
    is_ok <- x > 0 & x < 1
    if (any(!is_ok)) {
        obj <- deparse1(substitute(x))
        msg <- paste0(
            "All elements of ",
            if (val) "element `" else "argument `",
            obj,
            "` must be in (0, 1)."
        )
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}
