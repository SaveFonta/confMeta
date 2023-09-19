#' @title \emph{confMeta} objects
#' @description \code{confMeta} allows the user to create S3 objects of class
#'     \code{confMeta}.
#' @param estimates A vector containing the normalized individual effect
#'     estimates. Must be of the same length as \code{SEs} and coercible to
#'     type 'double'.
#' @param SEs The standard errors of the normalized individual effect estimates.
#'     Must be of the same length as \code{estimates} and coercible to
#'     type 'double'.
#' @param study_names Either \code{NULL} (the default) or a vector that can
#'     be coerced to type 'character'. Must be of the same length as
#'     arguments \code{estimates} and \code{SEs}. The vector is used to
#'     differentiate the individual studies. For more information, check out the
#'     Details section.
#' @param conf_level The confidence level. Must be a numeric vector of length
#'     one with it's value being in (0, 1).
#' @param fun A function that combines individual effect estimates and the
#'     corresponding standard errors into a combined p-value. The function must
#'     have arguments named 'estimates', 'SEs', and 'mu'.
#'     Additional arguments are also allowed, but these must either have a
#'     default value or be passed via the \code{...} argument.
#'     For more information, see the Details section.
#' @param ... Additional arguments passed to \code{fun}. See the Details
#'     section.
#' @return An S3 object of class \code{confMeta}. The object contains the
#'     following elements:
#'     \description{
#'         \item{estimates}{The normalized individual effect estimates.}
#'         \item{SEs}{The standard errors of the normalized individual effect
#'             estimates.}
#'         \item{study_names}{The names of the individual studies.}
#'         \item{conf_level}{The confidence level.}
#'         \item{individual_cis}{The confidence intervals for the
#'             individual effects. The exact calculation of these intervals
#'             can be found in the Details section.}
#'         \item{p_fun}{A function with arguments named 'estimates', 'SEs',
#'             'conf_level', and 'mu'. This is the p-value function that is
#'             used to find the combined confidence intervals.}
#'         \item{joint_cis}{The combined confidence intervals. The exact
#'             calculation of these intervals can be found in the Details
#'             section.}
#'     }
#' @details
#'     \section{Function arguments}{
#'     The argument \code{study_names} is used to differentiate between the
#'     different individual estimates. If the argument is set to \code{NULL},
#'     the element 'study_names' in the return object will just be the a
#'     character vector with elements "Study n" where n is a number from 1
#'     to \code{length(elements)}. These names are only used in some of the
#'     \code{autoplot} methods.
#'
#'     The argument \code{fun} must have arguments 'estimates', 'SEs',
#'     and 'mu' but it can also have further arguments.
#'     However, these must either have a default value or
#'     need to be passed via the \code{...} argument. If there
#'     are additional arguments passed via \code{...},
#'     \code{confMeta} will internally create a new
#'     function that calls \code{fun} with the additional arguments fixed.
#'     Since this is the p-value function that is used to calculate the
#'     combined confidence interval(s), it should return a vector of class
#'     'numeric' with value(s) in the interval [0, 1].
#'     }
#'
#'     \section{Confidence intervals}{
#'     The confidence intervals returned by \code{confMeta} are calculated as
#'     follows:
#'
#'     The \code{individual_cis} are calculated as
#'     \deqn{x_{i} \pm \Phi^{-1}{\text{conf_level}} \cdot \sigma_{i}}
#'     where \eqn{x_{i}} corresponds to the elements of vector
#'     \code{estimates}, \eqn{\Phi^{-1}} is the quantile function of
#'     the standard normal distribution, conf_level is the confidence
#'     level passed as argument \code{conf_level}, and
#'     \eqn{\sigma_{i}}, are the standard errors passed in argument
#'     \code{SEs}.
#'
#'     The \code{joint_cis} are found by searching where the function returned
#'     in element 'p_fun' is equal to 1-\code{conf_level}.
#'     }
#'
#' @examples
#'     # Simulate effect estimates and standard errors
#'     n <- 5
#'     estimates <- rnorm(n)
#'     SEs <- rgamma(n, 5, 5)
#'
#'     # Construct a simple confMeta object using p_edgington as
#'     # the p-value function
#'     cm <- confMeta(
#'         estimates = estimates,
#'         SEs = SEs,
#'         fun = p_edgington
#'     )
#'
#'     # Plot the object
#'     autoplot(cm, type = "p")                   # p-value function plot
#'     autoplot(cm, type = "forest")              # forest plot
#'     autoplot(cm, type = c("p", "forest"))      # both
#'
#' @export
confMeta <- function(
    estimates,
    SEs,
    study_names = NULL,
    conf_level = 0.95,
    fun,
    ...
) {

    # If study names is NULL, construct default names
    if (is.null(study_names)) {
        study_names <- paste0("Study ", seq_along(estimates))
    }

    # coerce inputs into correct format
    estimates <- as.double(estimates)
    SEs <- as.double(SEs)
    conf_level <- as.double(conf_level)
    study_names <- as.character(study_names)

    # run input checks
    validate_inputs(
        estimates = estimates,
        SEs = SEs,
        study_names = study_names,
        conf_level = conf_level,
        fun = fun
    )

    # Make the p-value function
    p_fun <- make_p_fun(
        fun = fun,
        ... = ...
    )

    new_confMeta(
        estimates = estimates,
        SEs = SEs,
        study_names = study_names,
        conf_level = conf_level,
        p_fun = p_fun
    )
}

# Constructor function
new_confMeta <- function(
    estimates = double(),
    SEs = double(),
    study_names = character(),
    conf_level = double(1L),
    p_fun,
    ...
) {

    # Calculate individual CIs
    alpha <- 1 - conf_level
    se_term <- alpha / 2 * SEs # Note: At some point include a one-sided option?
    individual_cis <- matrix(
        c(
            estimates - se_term,
            estimates + se_term
        ),
        ncol = 2L
    )

    # Calculate the joint CIs

    # Return object
    structure(
        list(
            estimates = estimates,
            SEs = SEs,
            conf_level = conf_level,
            p_fun = p_fun,
            individual_cis = individual_cis,
            joint_cis = joint_cis
        ),
        class = "confMeta"
    )
}

# Validator function
validate_inputs <- function(
    estimates,
    SEs,
    study_names,
    conf_level,
    fun
) {

    # Check inputs

    # type checks
    check_type(x = estimates, "double")
    check_type(x = SEs, "double")
    check_type(x = study_names, "character")
    check_type(x = conf_level, "double")
    check_type(x = fun, "closure")

    # Check lengths
    check_equal_length(                    # estimates, SEs, study_names should
        estimates = estimates,             # have same length
        SEs = SEs,
        study_names = study_names
    )
    check_length_1(x = conf_level)         # conf_level must be of length 1

    # Check validity of values
    check_all_finite(x = estimates)        # no NAs, NaNs etc in estimates
    check_all_finite(x = SEs)              # no NAs, NaNs etc in SEs
    check_all_finite(x = conf_level)       # no NAs, NaNs etc in conf_level
    check_prob(x = conf_level)             # conf_level must be between 0 & 1
    check_fun_args(x = fun)                # function must have correct args

    # Check that p_fun has an argument named mu

    # Return NULL if all checks passed
    invisible(NULL)
}

validate_confMeta <- function(confMeta) {

}


# ==============================================================================
# Make the p-value function  ===================================================
# ==============================================================================

make_p_fun <- function(fun, ...) {
    fa <- methods::formalArgs(fun)
    if (length(fa) > 4L && "mu" %in% fa) {
        # if the function is already only a function with exactly one argument
        # `mu`, return it as is
        fun
    } else {
        # otherwise construct a function that calls `fun` with the ellipsis
        # arguments fixed
        function(estimates, SEs, level, mu) {
            arglist <- append(list(...), alist(mu))
            do.call("fun", arglist)
        }
    }
}

################################################################################
# Checking the type of a variable                                              #
################################################################################


# ==============================================================================
# Argument checks ==============================================================
# ==============================================================================

################################################################################
# Checking the type of a variable                                              #
################################################################################

check_type <- function(x, type) {
    if (!(typeof(x) == type)) {
        obj <- deparse1(substitute(x))
        msg <- paste0("Argument `", obj, "` must be of type '", type, "'.")
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

################################################################################
# Checking the type of a variable                                              #
################################################################################

check_all_finite <- function(x) {
    if (!all(is.finite(x))) {
        obj <- deparse1(substitute(x))
        msg <- paste0(
            if (length(x) > 1L) {
                "All elements of argument `"
            } else {
                "Argument `"
            },
            obj,
            "` must be finite.",
        )
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

################################################################################
# Checking the type of a variable                                              #
################################################################################

check_type <- function(x, type) {
    if (!(typeof(x) == type)) {
        obj <- deparse1(substitute(x))
        msg <- paste0("Argument `", obj, "` must be of type '", type, "'.")
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

################################################################################
# Checking the type of a variable                                              #
################################################################################

check_class <- function(x, class) {
    if (!(class(x) == class)) {
        obj <- deparse1(substitute(x))
        msg <- paste0("Argument `", obj, "` must be of class '", class, "'.")
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

################################################################################
# Check that the function has correct arguments                                #
################################################################################

check_fun_args <- function(fun, ...) {

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
    # dotargs <- list(...)           # ellipsis args
    # return(dotargs)
    # nms_ellipsis <- names(dotargs) # names of the ellipsis args

    # Check whether the function has all the required
    ok_required <- has_all(
        x = must_have_args,
        superset = fa
    )

    # Check the non-required arguments whether they have a default
    # or are passed via ...
    non_standard <- fa[!(fa %in% must_have_args)]
    no_default <- vapply(
        f[non_standard],
        is_missing,
        logical(1L),
        USE.NAMES = TRUE
    )
    if (any(no_default)) {
        ell <- list(...)
        ell_args <- names(ell)

    }



    # Throw messages if there are unused arguments passed
    ok_unused <- has_unused(
        fun_args = fa,
        passed = c(must_have_args, nms_ellipsis)
    )

    if (length(dotargs) == 0) {
        # if there are no ellipsis args, check whether the function has
        # appropriate arguments
        ok <- all(fa %in% must_have_args) && all(must_have_args %in% fa)
        if (!ok) {
            msg <- paste0(
                "Argument `fun` must have arguments named ",
                paste0(
                    "'",
                    paste0(must_have_args, collapse = "', '"),
                    "'"
                ),
                "."
            )
            stop(msg)
        }

    }

    invisible(NULL)
}

dotargs <- check_fun_args(fun)

is_missing <- function(x) {
    if (!nzchar(x) && is.name(x)) {
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

check_equal_length <- function(...) {
    arguments <- list(...)
    l <- vapply(arguments, length, integer(1L), USE.NAMES = TRUE)
    ind <- l == l[1L]
    ok <- all(ind)
    if (ok) {
        invisible(NULL)
    } else {
        nms <- names(l)
        msg <- paste0(
            "Arguments ",
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

check_length_1 <- function(x) {
    if (length(x != 1L)) {
        obj <- deparse1(substitute(x))
        msg <- paste0("Argument `", obj, "` must be of length 1.")
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}

################################################################################
# Check for value between 0 and 1                                              #
################################################################################

check_prob <- function(x) {
    ok <- x > 0 & x < 1
    if (all(ok)) {
        obj <- deparse1(substitute(x))
        msg <- paste0("All elements of argument `", obj, "` must be in (0, 1).")
        stop(msg, call. = FALSE)
    }
    invisible(NULL)
}
