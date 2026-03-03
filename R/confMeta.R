#' @title Creating \code{confMeta} objects
#' @description Function to create objects of class \code{confMeta}. This is the
#'     main class within the package. For an overview of available methods
#'     run \code{methods(class = "confMeta")}.
#'     
#' @param estimates A vector containing the individual effect estimates. These should be on 
#' a scale where they are approximately normally distributed (e.g., log odds ratios, log risk 
#' ratios, mean differences).  Must be of the same length as \code{SEs} and coercible to 
#' type \code{double}.
#' @param SEs The standard errors of the individual effect estimates.
#'     Must be of the same length as \code{estimates} and coercible to
#'     type \code{double}.
#' @param study_names Optional vector of study identifiers. If \code{NULL} (the default) 
#'     names ("Study 1", "Study 2", ...) are used. Otherwise a vector that can
#'     be coerced to type \code{character}. Must be of the same length as
#'     arguments \code{estimates} and \code{SEs}. 
#' @param conf_level The confidence level (numeric scalar between 0 and 1).
#' @param fun A function that combines estimates and SEs into a combined \emph{p}-value. 
#'     The function must accept arguments \code{estimates}, \code{SEs}, and \code{mu}. 
#'     Any additional arguments that \code{fun} accepts can be passed via \code{...}. 
#'     See \code{\link{p_value_functions}} for supported methods.
#' @param fun_name A character vector of length 1. Label for \code{fun} used in
#'     \code{\link{autoplot.confMeta}}. Defaults to the name of the 
#'     function passed to \code{fun}.
#' @param w Numeric vector of weights for the studies. Must be of 
#'     the same length as \code{estimates} and \code{SEs}, finite, and non-negative. 
#'     If \code{NULL} (the default), equal weights are assumed. Used by weighted 
#'     \emph{p}-value combination functions.
#' @param MH Logical. If \code{TRUE}, the fixed-effect comparison method will use 
#'     Mantel-Haenszel pooling instead of the standard inverse-variance method. 
#'     Requires \code{table_2x2} and \code{measure} to be provided. Default is \code{FALSE}.
#' @param table_2x2 A data frame containing the 2x2 contingency table data for 
#'     Mantel-Haenszel pooling. Required if \code{MH = TRUE}. Must contain columns: 
#'     \code{ai} (events in experimental), \code{bi} (non-events in experimental), 
#'     \code{ci} (events in control), \code{di} (non-events in control), 
#'     \code{n1i} (total experimental), and \code{n2i} (total control).
#' @param measure A character string indicating the effect measure to be used 
#'     for Mantel-Haenszel pooling (e.g., "OR", "RR", "RD"). Required if 
#'     \code{MH = TRUE}.
#' @param method.tau Character string indicating which between-study variance 
#'     estimator to use for random-effects and Hartung-Knapp meta-analysis 
#'     (e.g., \code{"PM"}, \code{"REML"}, \code{"DL"}). Defaults to \code{"PM"}. 
#'     See \code{meta::metagen()} for all available choices.
#' @param adhoc.hakn.ci Character string indicating the variance correction 
#'     method for the Hartung-Knapp confidence intervals (e.g., \code{"IQWiG6"}, 
#'     \code{"HK"}). Defaults to \code{"IQWiG6"}. See \code{meta::metagen()} for choices.
#' @param ... Additional arguments passed to \code{fun}. 
#'     
#' @return
#' An S3 object of class \code{confMeta} containing the following elements:
#' \itemize{
#'  \item \code{estimates}, \code{SEs}, \code{w}, \code{study_names}, \code{conf_level}: Values 
#'  used to build the object.
#'  \item \code{individual_cis}: Matrix of Wald-type confidence intervals for each study.
#'   \item \code{p_fun}: The \emph{p}-value function used for combined inference.
#'   \item \code{joint_cis}: Combined confidence interval(s). Calculated
#'     by finding the values where the \emph{p}-value function is larger
#'     than the significant level (\code{1 - conf_level}).
#'   \item \code{gamma}: The local minima within the range of the individual effect
#'     estimates. Column \code{x} refers to the mean \eqn{\mu} and column \code{y}
#'     contains the corresponding \emph{p}-value.
#'   \item \code{p_max}: The local maxima of the \emph{p}-value function. Column \code{x}
#'     refers to the mean \eqn{\mu} and column \code{y} contains the corresponding
#'     \emph{p}-value.
#'   \item \code{p_0}: The value of the \emph{p}-value at \eqn{\mu = 0}.
#'   \item \code{comparison_cis}: Combined confidence intervals calculated with other
#'     methods. These can be used for comparison purposes. Currently, these
#'     other methods are fixed effect (IV or Mantel-Haenszel), random effects (IV), 
#'     Hartung & Knapp, and Henmi & Copas.
#'   \item \code{comparison_p_0}: The same as in element \code{p_0} but for the comparison
#'     methods (fixed effect, random effects, Hartung & Knapp, Henmi & Copas).
#'   \item \code{heterogeneity}: A data frame containing heterogeneity statistics 
#'     (Cochran's Q, \emph{p}-value for Q, \eqn{I^2}, \eqn{\tau^2}) calculated 
#'     using the REML estimator. 
#'   \item \code{table_2x2}: The 2x2 table data frame (if provided), otherwise \code{NULL}.
#' }
#'     
#' @section Confidence intervals:
#' Individual confidence intervals are calculated as:
#' \deqn{x_i \pm \Phi^{-1}(1 - \alpha/2) \cdot \sigma_i}
#' where \eqn{x_i} are the \code{estimates}, \eqn{\sigma_i} the \code{SEs}, and 
#' \eqn{\alpha = 1 - \mathrm{conf\_level}}.
#' 
#' Combined confidence intervals are found by inverting the \emph{p}-value function, 
#' identifying all \eqn{\mu} where the \emph{p}-value exceeds the significance level (1 - \code{conf_level}).
#' 
#' \strong{Note:} Due to necessity in further analysis, the method for estimating the between-study variance \eqn{\tau^2} in the Hartung-Knapp method for random effects meta-analysis uses the Paule-Mandel estimator and not the REML estimator anymore.
#' Also, the HK method uses \code{adhoc.hakn.ci = "IQWiG6"} by default, i.e., it uses variance correction if the HK confidence interval
#' is narrower than the CI from the classic random effects model with the DerSimonian-Laird estimator (IQWiG, 2022).
#' 
#' @examples
#'# Simulate effect estimates and standard errors
#'set.seed(42)
#'n <- 5
#'estimates <- rnorm(n)
#'SEs <- rgamma(n, 5, 5)
#'conf_level <- 0.95
#'
#'# Construct a simple confMeta object using p_edgington as
#'# the p-value function
#'cm <- confMeta(
#'   estimates = estimates,
#'   SEs = SEs,
#'   conf_level = conf_level,
#'   fun = p_edgington,
#'   fun_name = "Edgington",
#'   input_p = "greater"
#')
#'     
#'cm2 <- confMeta(
#'   estimates = estimates,
#'   SEs = SEs,
#'   conf_level = conf_level,
#'   fun = p_edgington_w,
#'   w = 1/SEs,
#'   fun_name = "Edgington (1/SE)",
#'   input_p = "greater"
#')
#' 
#'# print the objects
#'cm
#'cm2
#' 
#'# Plot the objects
#'autoplot(cm, cm2, type = "p")                   # p-value function plot
#'autoplot(cm, cm2, type = "forest")              # forest plot
#'autoplot(cm, cm2, type = c("p", "forest"))      # both
#'
#' @seealso \code{\link{p_value_functions}}, \code{\link{autoplot.confMeta}}
#' @export
confMeta <- function(
    estimates,
    SEs,
    study_names = NULL,
    conf_level = 0.95,
    fun,
    fun_name = NULL,
    w = NULL,    
    MH = FALSE, 
    table_2x2 = NULL,
    measure = NULL,
    method.tau = "PM", 
    adhoc.hakn.ci = "IQWiG6",
    ...
) {
  
  # If study names is NULL, construct default names
  if (is.null(study_names)) {
    study_names <- paste0("Study ", seq_along(estimates))
  }
  # If the function name is not given, add a default one
  if (is.null(fun_name)) {
    fun_name <- deparse1(substitute(fun))
  }
  
  # Catch the ... arguments
  ell <- list(...)
  ell <- remove_unused(fun = fun, ell = ell)
  
  
  # coerce inputs into correct format
  if (inherits(estimates, "numeric")) estimates <- as.double(estimates)
  if (inherits(SEs, "numeric")) SEs <- as.double(SEs)
  if (!is.null(w)) w <- as.double(w) #[MOD]
  if (inherits(conf_level, "numeric")) conf_level <- as.double(conf_level)
  study_names <- as.character(study_names)
  
  
  # run input checks
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
  
  
  # Add input checks for MH
  if (MH) {
    if (is.null(table_2x2)) {
      stop("When MH = TRUE, 'table_2x2' must be provided.")
    }
    if (is.null(measure)) {
      stop("When MH = TRUE, 'measure' must be provided.")
    }
    # Force data.frame, otw matrix can break this
    if (!is.data.frame(table_2x2)) {
      table_2x2 <- as.data.frame(table_2x2)
    }
    
    
    # Check table_2x2 structure
    required_cols <- c("ai", "bi", "ci", "di", "n1i", "n2i")
    missing_cols <- setdiff(required_cols, names(table_2x2))
    if (length(missing_cols) > 0) {
      stop("table_2x2 must contain columns: ", 
           paste(missing_cols, collapse = ", "))
    }
  }
  
  
  # Make the p-value function
  p_fun <- make_p_fun(
    fun = fun,
    ell = ell
  )
  
  new_confMeta(
    estimates = estimates,
    SEs = SEs,
    w = w,  # [MOD] can be NULL, new_confMeta will take care
    study_names = study_names,
    conf_level = conf_level,
    p_fun = p_fun,
    fun_name = fun_name,
    MH = MH, 
    table_2x2 = table_2x2,
    measure = measure,
    method.tau = method.tau,
    adhoc.hakn.ci = adhoc.hakn.ci
  )
}

# Constructor function
#' @importFrom stats qnorm
#' @importFrom meta metagen metabin
#' @noRd

new_confMeta <- function(
    estimates = double(),
    SEs = double(),
    w = NULL,   # [MOD] 
    study_names = character(),
    conf_level = double(1L),
    p_fun,
    fun_name,
    MH = FALSE, 
    table_2x2 = NULL,
    measure = NULL,
    method.tau = "PM",
    adhoc.hakn.ci = "IQWiG6"
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
    w = w,  # [MOD] if not null, will be used
    conf_level = conf_level,
    p_fun = p_fun
  )
  
  # Calculate joint CIs with the comparison methods
  method <- c("fe", "re", "hk", "hc")
  comparison <- get_stats_others(
    method = method,
    estimates = estimates,
    SEs = SEs,
    conf_level = conf_level,
    method.tau = method.tau,      
    adhoc.hakn.ci = adhoc.hakn.ci 
  )
  
  # overwrite the FE method if MH = TRUE 
  if (MH == TRUE) comparison <- overwrite_FE (comparison = comparison, table_2x2 = table_2x2, measure = measure, conf_level = conf_level) 
  
  metagen_obj <- metagen_wrap(estimates, SEs, conf_level = conf_level)
  
  
  heterogeneity <- data.frame(
    Q = metagen_obj$Q,
    p_Q = metagen_obj$pval.Q,
    I2 = metagen_obj$I2,
    Tau = metagen_obj$tau,
    I2_lower = metagen_obj$lower.I2,
    I2_upper = metagen_obj$upper.I2,
    stringsAsFactors = FALSE
  )
  
  
  
  # Return object
  structure(
    list(
      estimates = estimates,
      SEs = SEs,
      w = w,   # [MOD] 
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
      comparison_p_0 = comparison$p_0,
      heterogeneity = heterogeneity,
      table_2x2 = table_2x2
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
    "w",   
    "heterogeneity",
    "table_2x2"  
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
      
      if (!is.null(table_2x2)) {
        check_class(x = table_2x2, class = "data.frame", val = TRUE)
      }
      
      invisible(NULL)
    }
  )
  
  invisible(NULL)
}

# ==============================================================================
# Calculate the CIs using the other methods ====================================
# ==============================================================================

#' @importFrom meta metagen
#' @noRd
get_obj_re <- function(estimates, SEs, conf_level, method.tau) {
    meta::metagen(
        TE = estimates, seTE = SEs, sm = "MD",
        level = conf_level, method.tau = method.tau,
        random = TRUE, common = FALSE
    )
}
#???
#' @importFrom meta metagen
#' @noRd
get_obj_fe <- function(estimates, SEs, conf_level) {
    meta::metagen(
        TE = estimates, seTE = SEs, sm = "MD",
        level = conf_level, method.tau = "REML", #note that here method.tau doesn't change anything 
        random = FALSE, common = TRUE
    )
}

#' @importFrom meta metagen
#' @noRd
get_obj_hk <- function(estimates, SEs, conf_level, method.tau, adhoc.hakn.ci) {
    meta::metagen(
        TE = estimates, seTE = SEs, sm = "MD",
        level = conf_level, method.tau = method.tau, method.random.ci = "HK", adhoc.hakn.ci = adhoc.hakn.ci, #[MOD]--> "hakn = TRUE" is deprecated
        common = FALSE, random = TRUE
    )
} # IMPORTANT --> reading the documentation of metagen we have that method.tau = gs("method.tau"), this mean that by default it will
# try to find a method.tau in the env, that can be actually modified using settings.meta(method.tau = "PM"), even though
# it can be useful, this would mean that at the beginning of every confMeta session, we should set the method, and it is not
# always practical and straightforward to new users of the package. So I decided to force it to PM for the moment.
# by default when loading the package, the general setting is REML (can see this using meta::settings.meta()$method.tau)

#I also added the ad.hoc correction

#' @importFrom metafor hc rma
#' @noRd
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
#' @noRd
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
get_stats_others <- function(method, estimates, SEs, conf_level, method.tau, adhoc.hakn.ci) {
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
                    conf_level = conf_level,
                    method.tau = method.tau
                ),
                "hk" = get_obj_hk(
                    estimates = estimates,
                    SEs = SEs,
                    conf_level = conf_level,
                    method.tau = method.tau,
                    adhoc.hakn.ci = adhoc.hakn.ci
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


# function to overwrite the FE wit MH if 2x2 table provided
overwrite_FE <- function(comparison, table_2x2, measure, conf_level = 0.95) {
  tryCatch({  
    # alternative implementation with rma.mh, problem --> doesnt handle 0 cells properly
    #MH_fe <- metafor::rma.mh(
    #  ai = table_2x2$ai,
     # bi = table_2x2$bi,
    #  ci = table_2x2$ci,
    #  di = table_2x2$di,
    #  n1i = table_2x2$n1i,
    #  n2i = table_2x2$n2i,
    #  measure = measure,
    #  level = conf_level * 100  # rma.mh wants percentage
    #)
    #comparison$CI["Fixed effect", "lower"] <- MH_fe$ci.lb
    #comparison$CI["Fixed effect", "upper"] <- MH_fe$ci.ub
    #comparison$p_0["Fixed effect", "y"] <- MH_fe$pval 
    
    MH_fe <- meta::metabin(
      event.e =  table_2x2$ai,           # events in experimental group
      n.e = table_2x2$n1i,              # Number of obs in experimental group 
      event.c = table_2x2$ci,          # events in control group   
      n.c = table_2x2$n2i,              #  obs in control group 
      sm = measure,              # Summary Measure
      method = "MH",
      allstudies = TRUE, # include studies with zero or all events in both groups 
      level = conf_level
      )
    
    # NOTE FOR THE FUTURE --> In a revision, I was scared that I needed to add backtransf = TRUE, but backtranf is just related to
    #printouts and plots, so CIs are reported in the right scale
    
    comparison$CI["Fixed effect", "lower"] <- MH_fe$lower.common
    comparison$CI["Fixed effect", "upper"] <- MH_fe$upper.common
    comparison$p_0["Fixed effect", "y"] <- MH_fe$pval.common 
 
    comparison
  }, error = function(e) {
    warning("Mantel-Haenszel pooling failed: ", e$message, 
            "\nFalling back to inverse-variance FE method.")
    comparison  # Return unchanged
  })
}




################################################################################
# Compute general metagen object                                               #
################################################################################

# with this function we compute a general metagen object. If one need to add other statistics computed by metagen (since it computes
#everything you can think of in meta analysis), just extract from this (NB: here you can easily change the default options used in metagen to estimate tau in different ways)

# NOTE: I am aware that "metagen_wrap" and "get_obj_re" do the same call to "metagen"
# But I wanted them to keep them separate to have more control on the estimator used in each one of them 
# if in the future we want to change them 
#???
metagen_wrap <- function(estimates, SEs, conf_level = 0.95) {
  meta::metagen(
    TE = estimates, seTE = SEs,
    level = conf_level,
    method.tau = "REML",
    random = TRUE, common = FALSE
  )
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




