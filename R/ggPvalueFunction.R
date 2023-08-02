#' Plot the p-value function(s) for a given set of studies
#'
#' @template thetahat
#' @template se
#' @template level
#' @param distr Character vector of length 1 or 2. Valid options are
#' either \code{"f"} or \code{"chisq"} or both. Denotes the distribution
#' used in cases where \code{pValueFUN = "hMean"} but ignored otherwise.
#' @param pValueFUN Character vector of length 1, 2 or 3. Valid options consist
#' of any combination of \code{"hMean"}, \code{"k-Trials"}, or \code{"Pearson"}.
#' \code{"hMean"} will add a line using \code{\link[hMean]{hMeanChiSqMu}} as the
#' p-value function, \code{"k-Trials"} adds a line using
#' \code{\link[hMean]{kTRMu}} as the p-value function, and \code{"Pearson} uses
#' \code{\link[hMean]{pPearsonMu}} as the p-value function.
#' @param heterogeneity The heterogeneity model used to adjust the standard
#' errors \code{se}. Valid options are any combination of \code{"none"},
#' \code{"additive"}, or \code{"multiplicative"}. See also
#' \code{\link[hMean]{hMeanChiSqMu}}.
#' @param pValueFUN Any combination of \code{"hMean"}, \code{"k-Trials"}, and
#' \code{"Pearson"}. This argument determines which methods are plotted in the
#' resulting image.
#' @template pValueFUN_args
#' @param xlim Numeric vector of length 2 denoting the range of values
#' on the x-axis, i.e. the effect size.
#'
#' @return An object of class \code{ggplot}. It contains the plot of the p-value
#' functions specified by the input arguments.
#'
#' @examples
#' thetahat <- c(-0.49, -0.17, -0.52, -0.48, -0.26, -0.36, -0.47, -0.3, -0.15,
#'               -0.28)
#' se <- c(0.346945150708765, 0.214289651908355, 0.239800324754587,
#'         0.372455823554997, 0.25000459389308, 0.295923805016299,
#'         0.219391786477601, 0.14796190250815, 0.270413132170067,
#'         0.500009187786161)
#' pValueFUN <- c("hMean", "k-Trials", "Pearson", "Edgington")
#' distr <- c("chisq")
#' heterogeneity <- "none"
#' ggPvalueFunction(thetahat = thetahat,
#'                  se = se,
#'                  xlim = c(-0.8, 0.1),
#'                  distr = distr,
#'                  heterogeneity = heterogeneity,
#'                  pValueFUN = pValueFUN)
#' @export
ggPvalueFunction <- function(
    thetahat,
    se,
    level = 0.95,
    distr = c("chisq", "f"),
    heterogeneity = c("none", "additive", "multiplicative"),
    pValueFUN = c("hMean", "k-Trials", "Pearson", "Edgington", "Fisher"),
    pValueFUN_args,
    xlim = c(min(thetahat - 2 * se), max(thetahat + 2 * se))
) {

    # get the p-value function(s)
    pValueFUN <- match.arg(pValueFUN, several.ok = TRUE)

    # get the desired heterogeneity
    heterogeneity <- match.arg(heterogeneity, several.ok = TRUE)

    # get the distribution
    distr <- match.arg(distr, several.ok = TRUE)

    # Check the other inputs
    stopifnot(
        is.numeric(level),
        level > 0 && level < 1,
        is.numeric(thetahat),
        length(thetahat) > 0L,
        is.numeric(se),
        length(se) == length(thetahat) || length(se) == 1L,
        is.numeric(xlim),
        length(xlim) == 2L,
        xlim[1] < xlim[2]
    )

    # Construct the grid to loop over
    grid <- make_grid(
        pValueFUN = pValueFUN,
        heterogeneity = heterogeneity,
        distr = distr
    )

    # Set some constants that are equal for all grid rows
    const <- list(
        thetahat = thetahat,
        se = se,
        alpha = 1 - level,
        phi = estimatePhi(thetahat = thetahat, se = se),
        tau2 = estimateTau2(
            thetahat = thetahat,
            se = se,
            control = list(stepadj = 0.5, maxiter = 1000, threshold = 1e-6)
        ),
        eps = 0.0025, # for plotting error bars
        eb_height = 0.025,
        muSeq = seq(xlim[1], xlim[2], length.out = 1e4)
    )

    # Calculate the p-values and CIs
    data <- lapply(seq_len(nrow(grid)), function(x) {
        grid_row <- grid[x, ]
        p_call <- make_p_call(grid_row = grid_row, const = const)
        CI_call <- make_CI_call(p_call = p_call, level = level)
        pval <- eval(p_call)
        CIs <- eval(CI_call)

        # # make a data frame that contains the necessary information
        # # gamma_min (for point on the minimum), p-values, etc.
        # idx <- which.min(CIs$gamma[, 2])
        # gamma_exists <- length(idx) != 0L
        # gamma_min <- if (!gamma_exists) NA_real_ else CIs$gamma[idx, 2]
        # x_gamma_min <- if (!gamma_exists) NA_real_ else CIs$gamma[idx, 1]

        # Instead of gamma_min, calculate the value where the p-value function
        # is 0
        y0_call <- p_call
        y0_call[["mu"]] <- 0
        y0 <- eval(y0_call)

        df1 <- data.frame(
            x = const$muSeq,
            y = pval,
            heterogeneity = grid_row$heterogeneity,
            p_val_fun = grid_row$fun_name,
            distr = grid_row$distr,
            y0 = y0,
            group = factor(grid_row$name, levels = grid$name),
            stringsAsFactors = FALSE,
            row.names = NULL
        )
        # handle error bars (display errorbars with a little jitter depending on
        # heterogeneity)
        factor <- switch(
          grid_row$heterogeneity,
          "none" = 0,
          "additive" = 1,
          "multiplicative" = -1
        )
        # make a second data frame for the display of confidence intervals
        df2 <- data.frame(
            xmin = CIs$CI[, 1],
            xmax = CIs$CI[, 2],
            heterogeneity = grid_row$heterogeneity,
            p_val_fun = grid_row$fun_name,
            distr = grid_row$distr,
            y = rep(1 - level, nrow(CIs$CI)) + factor * const$eps,
            group = factor(grid_row$name, levels = grid$name),
            stringsAsFactors = FALSE,
            row.names = NULL
        )
        df2$ymax <- df2$y + const$eb_height
        df2$ymin <- df2$y - const$eb_height
        list(df1, df2)
    })
    # Extract the data frame for the lines with p-value functions
    # as well as the data frame for the error bars
    plot_data <- lapply(
        list(lines = 1L, errorbars = 2L),
        function(z, data) do.call("rbind", lapply(data, "[[", i = z)),
        data = data
    )
    lines <- plot_data[["lines"]]
    errorbars <- plot_data[["errorbars"]]

    # Define function to convert breaks from primary y-axis to
    # breaks for secondary y-axis
    trans <- function(x) rev(abs(x - 1) * 100)
    # Define breaks for the primary y-axis
    breaks_y1 <- sort(c(1 - level, pretty(c(lines$y, 1))))
    # Compute breaks for the secondary y-axis
    breaks_y2 <- trans(breaks_y1)
    # Set transparency
    transparency <- 1

    p <- ggplot2::ggplot(
        data = lines,
        ggplot2::aes(x = x, y = y, color = group)
    ) +
    ggplot2::geom_hline(yintercept = 1 - level, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = thetahat, linetype = "dashed") +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid") +
    ggplot2::geom_line(alpha = transparency) +
    ggplot2::geom_point(
        data = lines[!is.na(lines$y0), ],
        ggplot2::aes(x = 0, y = y0, color = group),
        alpha = transparency
    ) +
    ggplot2::scale_y_continuous(
        name = "p-value",
        breaks = breaks_y1,
        limits = c(0, 1),
        sec.axis = ggplot2::sec_axis(
            trans = trans,
            name = "Confidence level [%]",
            breaks = breaks_y2
        )
    ) +
    # Draw intervals on x-axis
    ggplot2::geom_segment(
        data = errorbars[!is.na(errorbars$xmin), ],
        ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y)
    ) +
    ggplot2::geom_segment(
        data = errorbars[!is.na(errorbars$xmin), ],
        ggplot2::aes(x = xmin, xend = xmin, y = ymin, yend = ymax)
    ) +
    ggplot2::geom_segment(
        data = errorbars[!is.na(errorbars$xmin), ],
        ggplot2::aes(x = xmax, xend = xmax, y = ymin, yend = ymax)
    ) +
    # Set x-axis window, labels
    ggplot2::xlim(xlim) +
    ggplot2::labs(
        x = bquote(mu),
        color = "Configuration"
    ) +
    # Set theme
    ggplot2::theme_minimal() +
    ggplot2::theme(
        axis.title.y.right = ggplot2::element_text(angle = 90),
        legend.position = "bottom"
    )

    # make the plot
    print(p)

    # return
    invisible(
        list(plot = p, p_0 = unique(lines[c("group", "y0")]))
    )
}

## This function constructs calls
make_p_call <- function(grid_row, const) {
    # Check inputs
    if (!inherits(grid_row, "data.frame"))
        stop("Argument 'grid_row' must be a data.frame.")
    grid_names <- c("fun_name", "heterogeneity", "distr")
    if (!all(grid_names %in% names(grid_row)))
        stop(
            paste0(
                "Argument 'grid_row' must have names: ",
                paste0(grid_names, collapse = ", ")
            )
        )
    if (nrow(grid_row) != 1L)
        stop("Argument 'grid_row' must have exactly 1 row.")
    # Construct the call to the p-value function
    currentFUN <- get(grid_row$fun_name, pos = "package:hMean")
    args <- list(
        thetahat = const$thetahat,
        se = const$se,
        heterogeneity = grid_row$heterogeneity,
        phi = if (grid_row$heterogeneity == "multiplicative")
            const$phi
        else NULL,
        tau2 = if (grid_row$heterogeneity == "additive")
            const$tau2
        else NULL,
        mu = const$muSeq
    )
    ## If pvalueFUN = hMean, add the distribution
    if (grid_row$fun_name == "hMeanChiSqMu") {
        args <- append(
            args,
            list(distr = grid_row$distr),
        )
    }
    ## Add the default arguments
    def_args <- formals(currentFUN)
    def_args <- def_args[!names(def_args) %in% names(args)]
    ## Force evaluation of default args
    def_args <- with(args, lapply(def_args, function(x) eval(x)))
    ## Add default arguments to actual arguments
    args <- append(args, def_args)
    as.call(append(list(currentFUN), args))
}

make_CI_call <- function(p_call, level) {

    # convert p_call into list
    p_call_list <- as.list(p_call)

    # remove function, thetahat, se, mu as they are passed anyway
    remove <- names(p_call_list) %in% c("", "thetahat", "se", "mu")
    pValueFUN_args <- p_call_list[!remove]

    # add the function as pValueFUN
    pValueFUN <- p_call_list[[1L]]

    # Assemble the other args
    as.call(
      list(
        hMeanChiSqCI,
        thetahat = p_call_list$thetahat,
        se = p_call_list$se,
        level = level,
        pValueFUN = pValueFUN,
        pValueFUN_args = pValueFUN_args
      )
    )
}
