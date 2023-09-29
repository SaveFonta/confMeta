autoplot.confMeta <- function(
    estimates,
    SEs,
    conf_level = 0.95,
    xlim = c(min(thetahat - 2 * se), max(thetahat + 2 * se)),
    drapery = TRUE
) {

    # Check the other inputs
    stopifnot(
        is.numeric(conf_level),
        conf_level > 0 && conf_level < 1,
        is.numeric(estimates),
        length(estimates) > 0L,
        is.numeric(SEs),
        length(SEs) == length(estimates) || length(SEs) == 1L,
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
        estimates = estimates,
        SEs = SEs,
        alpha = 1 - conf_level,
        phi = estimate_phi(estimates = estimates, SEs = SEs),
        tau2 = estimate_tau2(
            estimates = estimates,
            SEs = SEs,
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
        CI_call <- make_CI_call(p_call = p_call, conf_level = conf_level)
        pval <- eval(p_call)
        CIs <- eval(CI_call)

        # hMeanChiSqCI(
        #     thetahat = CI_call[["thetahat"]],
        #     se = CI_call[["se"]],
        #     level = CI_call[["level"]],
        #     pValueFUN = CI_call[["pValueFUN"]],
        #     pValueFUN_args = CI_call[["pValueFUN_args"]],
        #     alternative = "none",
        #     check_inputs = TRUE
        # )

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
            y = rep(1 - conf_level, nrow(CIs$CI)) + factor * const$eps,
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

    # Calculate the drapery lines
    if (drapery) {
        dp <- get_drapery_df(
            estimates = const$estimates,
            SEs = const$SEs,
            mu = const$muSeq
        )
    }

    # Define function to convert breaks from primary y-axis to
    # breaks for secondary y-axis
    trans <- function(x) abs(x - 1) * 100
    # Define breaks for the primary y-axis
    b_points <- c(1 - conf_level, pretty(c(lines$y, 1)))
    o <- order(b_points, decreasing = FALSE)
    breaks_y1 <- b_points[o]
    # Compute breaks for the secondary y-axis
    # breaks_y2 <- trans(b_points[o])
    # Set transparency
    transparency <- 1

    p <- ggplot2::ggplot(
        data = lines,
        ggplot2::aes(x = x, y = y, color = group)
    ) +
    ggplot2::geom_hline(yintercept = 1 - conf_level, linetype = "dashed")
    if (!drapery) {
        p <- p + ggplot2::geom_vline(
            xintercept = estimates,
            linetype = "dashed"
        )
    } else {
        p <- p + ggplot2::geom_line(
            data = dp,
            mapping = ggplot2::aes(x = x, y = y, group = study),
            linetype = "dashed",
            color = "lightgrey",
            show.legend = FALSE
        )
    }
    p <- p +
    ggplot2::geom_vline(xintercept = 0, linetype = "solid") +
    ggplot2::geom_hline(yintercept = 0, linetype = "solid") +
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
        expand = c(0, 0),
        sec.axis = ggplot2::sec_axis(
            trans = trans,
            name = "Confidence level [%]",
            breaks = trans(breaks_y1)
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

    # Return also the value at 0
    p_0 <- unique(lines[c("group", "y0")])
    names(p_0) <- c("function", "p_0")
    rownames(p_0) <- NULL

    # return
    list(
        plot = p,
        p_0 = p_0
    )
}

# Calculate the drapery lines
get_drapery_df <- function(estimates, SEs, mu) {
    # get lenghts
    l_t <- length(estimates)
    l_m <- length(mu)
    l_tot <- l_t * l_m
    # Initialize vectors
    x_dp <- rep(mu, times = l_t)
    y_dp <- study <- numeric(l_tot)
    # Indices to loop over
    idx <- seq_len(l_m)
    for (i in seq_along(estimates)) {
        y_dp[idx] <- 2 *
        (1 - stats::pnorm(abs(estimates[i] - mu) / SEs[i]))
        study[idx] <- rep(i, l_m)
        idx <- idx + l_m
    }
    data.frame(
        x = x_dp,
        y = y_dp,
        study = study,
        stringsAsFactors = FALSE
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
    currentFUN <- get(grid_row$fun_name, pos = "package:confMeta")
    args <- list(
        estimates = const$estimates,
        SEs = const$SEs,
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

make_CI_call <- function(p_call, conf_level) {

    # convert p_call into list
    p_call_list <- as.list(p_call)

    # remove function, thetahat, se, mu as they are passed anyway
    remove <- names(p_call_list) %in% c("", "estimates", "SEs", "mu")
    pValueFUN_args <- p_call_list[!remove]

    # add the function as pValueFUN
    pValueFUN <- p_call_list[[1L]]

    # Assemble the other args
    as.call(
      list(
        hMeanChiSqCI,
        estimates = p_call_list$estimates,
        SEs = p_call_list$SEs,
        conf_level = conf_level,
        pValueFUN = pValueFUN,
        pValueFUN_args = pValueFUN_args
      )
    )
}
