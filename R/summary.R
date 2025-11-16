#' @title Summarize one or more confMeta objects
#'
#' @description
#' Produces a compact tabular summary of one or more `confMeta` objects,
#' including the combined confidence intervals, \eqn{p_0}, and (if available)
#' results from the comparison methods.
#'
#' When multiple `confMeta` objects are supplied (e.g. with different
#' combination functions), their summaries are stacked in a single table.
#' The comparison methods are taken from the first object.
#'
#'
#'
#' @param object A `confMeta` object, as returned by [confMeta()].
#' @param ... Additional arguments (currently ignored).
#' 
#' 
#' @return
#' A data frame where each row corresponds to a method:
#'
#' - the main method of each `confMeta` object (named according to `fun_name`)
#' - comparison methods (e.g. fixed-effect, random-effects, Hartung–Knapp,
#'   Henmi–Copas), taken from the first object
#'
#' with the following columns:
#'
#' - `upper`, `lower` – bounds of the confidence interval  
#' - `AUCC`, `AUCC_ratio` – available only for the main methods (one row per
#'   `confMeta` object)  
#' - `p_0` – p-value at \eqn{\mu = 0}
#'
#'
#' @export



#' @export
summary.confMeta <- function(object, ...) {
  # collect all confMeta objects
  objs <- c(list(object), list(...))
  
  is_confMeta <- vapply(objs, inherits, logical(1), what = "confMeta")
  if (!all(is_confMeta)) {
    stop("All arguments to summary.confMeta must be objects of class 'confMeta'.")
  }
  
  #maybe one day I will add also the bayesmeta, now I don't have time 
  
  # extract main row from a single confMeta (if the method produces two CI, not supported, but shouldn't happen)
  extract_main_row <- function(obj) {
    main_upper <- main_lower <- main_p0 <- main_aucc <- main_aucc_ratio <- NA_real_
    main_name <- obj$fun_name
    
    # joint_cis: take first interval (if multiple)
    if (!is.null(obj$joint_cis)) {
      jc <- as.data.frame(obj$joint_cis)
      if (!all(c("lower", "upper") %in% names(jc))) {
        stop("`joint_cis` must have columns named 'lower' and 'upper'.")
      }
      main_lower <- jc$lower[1L]
      main_upper <- jc$upper[1L]
    }
    
    # p_0: 1x2 (x, y); take y as p-value
    if (!is.null(obj$p_0)) {
      p0 <- as.data.frame(obj$p_0)
      if ("y" %in% names(p0)) {
        main_p0 <- p0$y[1L]
      }
    }
    
    if (!is.null(obj$aucc))       main_aucc       <- obj$aucc
    if (!is.null(obj$aucc_ratio)) main_aucc_ratio <- obj$aucc_ratio
    
    data.frame(
      upper      = main_upper,
      lower      = main_lower,
      AUCC       = main_aucc,
      AUCC_ratio = main_aucc_ratio,
      p_0        = main_p0,
      row.names  = main_name
    )
  }
  
  # extract comparison rows from a single confMeta 
  extract_comparisons <- function(obj) {
    if (is.null(obj$comparison_cis)) {
      return(NULL)
    }
    
    comp_cis <- as.data.frame(obj$comparison_cis)
    if (!all(c("lower", "upper") %in% names(comp_cis))) {
      stop("`comparison_cis` must have columns named 'lower' and 'upper'.")
    }
    
    method_names <- rownames(comp_cis)
    if (is.null(method_names) || any(method_names == "")) {
      method_names <- paste0("method_", seq_len(nrow(comp_cis)))
      rownames(comp_cis) <- method_names
    }
    
    # p_0 for comparisons: take column y from comparison_p_0 if present
    p0_comp <- rep(NA_real_, length(method_names))
    if (!is.null(obj$comparison_p_0)) {
      comp_p0 <- as.data.frame(obj$comparison_p_0)
      if ("y" %in% names(comp_p0)) {
        p0_comp <- comp_p0[method_names, "y"]
      }
    }
    
    data.frame(
      upper      = comp_cis$upper,
      lower      = comp_cis$lower,
      AUCC       = NA_real_,
      AUCC_ratio = NA_real_,
      p_0        = p0_comp,
      row.names  = method_names
    )
  }
  
  # main rows: one per confMeta object 
  main_list <- lapply(objs, extract_main_row)
  main_df   <- do.call(rbind, main_list)
  
  # comparison rows: taken from the first object only
  comp_df <- extract_comparisons(objs[[1L]])
  
  # combine and return
  if (is.null(comp_df)) {
    out <- main_df
  } else {
    out <- rbind(main_df, comp_df)
  }
  
  out
}
