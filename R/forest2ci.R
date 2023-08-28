#' Extract the confidence interval from forest plots
#'
#' @param forest The output of a call to \code{\link[confMeta]{ForestPlot}}
#'
#' @return A list with confidence intervals from all the methods specified
#' in the call to \code{ForestPlot}.
#'
#' @examples
#' thetahat <- c(-0.49, -0.17, -0.52, -0.48, -0.26, -0.36, -0.47, -0.3, -0.15,
#'               -0.28)
#' se <- c(0.346945150708765, 0.214289651908355, 0.239800324754587,
#'         0.372455823554997, 0.25000459389308, 0.295923805016299,
#'         0.219391786477601, 0.14796190250815, 0.270413132170067,
#'         0.500009187786161)
#' pValueFUN <- c("hMean", "k-Trials", "Pearson", "Edgington", "Fisher")
#' distr <- c("chisq")
#' heterogeneity <- "none"
#' forest <- ForestPlot(
#'     thetahat = thetahat,
#'     se = se,
#'     distr = distr,
#'     pValueFUN = pValueFUN,
#'     heterogeneity = heterogeneity
#' )
#' forest2CI(forest)
#'
#' @export
forest2CI <- function(forest) {
    data <- forest[["plot"]][["plot_env"]]$new_method_cis[["CIs"]]
    s <- split(data, f = data$name)
    lapply(s, get_lower_upper)
}

################################################################################
# Helper function to extract the CI for one method                             #
################################################################################
get_lower_upper <- function(method_table) {
    if (all(is.na(method_table$id))) {
        out <- matrix(
            rep(NA_real_, 2L),
            ncol = 2L,
            dimnames = list(NULL, c("lower", "upper"))
        )
    } else {
        ss <- split(method_table, method_table$id)
        out <- t(
            vapply(
                ss,
                function(dat) with(dat, x[c(1L, length(x) / 2L + 1L)]),
                double(2L),
                USE.NAMES = FALSE
            )
        )
        colnames(out) <- c("lower", "upper")
        out
    }
}
