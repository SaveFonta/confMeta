#' @export
#' @noRd
print.confMeta <- function(x, ...) {
  
  k <- length(x$estimates)
  cat("  Meta-Analysis with p-value combination\n")
  cat("-----------------------------------------\n")
  cat("Number of studies: ", k, "\n")
  cat("Combination method:", x$fun_name, "\n\n")
  
  point_est <- round(x$p_max[1, "x"], 3)
  conf_pct  <- x$conf_level * 100
  
  lower <- round(x$joint_cis[, "lower"], 3)
  upper <- round(x$joint_cis[, "upper"], 3)
  
  cat(sprintf("Estimate: %.3f,\n%.0f%% Confidence Interval: [%.3f, %.3f]\n", 
              point_est, conf_pct, lower, upper), sep = "")
  
  cat("\n")
  invisible(x)
}
