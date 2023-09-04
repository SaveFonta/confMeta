# Constructor function
new_confMeta <- function(
    estimates = double(),
    SEs = double(),
    conf_level = double(1L),
    p_fun,
    ...
) {

    # Calculate individual CIs
    alpha <- 1 - conf_level
    se_term <- alpha / 2 * SEs # Maybe include a one-sided option?
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
    study_names = NULL,
    conf_level = 0.95,
    p_fun
) {

    # Check inputs

    # type checks
    check_type(x = estimates, "double")
    check_type(x = SEs, "double")
    check_type(x = conf_level, "double")
    check_type(x = p_fun, "closure")

    # Check lengths
    check_equal_length(                    # estimates and SEs of same length
        estimates = estimates,
        SEs = SEs
    )
    check_length_1(x = conf_level)         # conf_level must be of length 1

    # Check validity of values
    check_all_finite(x = estimates)        # no NAs, NaNs etc in estimates
    check_all_finite(x = SEs)              # no NAs, NaNs etc in SEs
    check_all_finite(x = conf_level)       # no NAs, NaNs etc in conf_level
    check_prob(x = conf_level)             # conf_level must be between 0 & 1

    # Return NULL if all checks passed
    invisible(NULL)
}

validate_confMeta <- function(confMeta) {

}


# Helper function
confMeta <- function(
    estimates,
    SEs,
    studyNames = NULL,
    conf_level = 0.95,
    p_fun,
    ...
) {

    # coerce inputs into correct format
    estimates <- as.double(estimates)
    SEs <- as.double(SEs)
    conf_level <- as.double(conf_level)

    # recycle to appropriate lengths

    # run input checks
    validate_inputs(
        estimates = estimates,
        SEs = SEs,
        study_names = study_names,
        conf_level = conf_level,
        p_fun = p_fun
    )

    validate_confMeta(
        new_confMeta(
            estimates = estimates,
            SEs = SEs,
            conf_level = conf_level,
            p_fun = p_fun,
            ... = ...
        )
    )
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
