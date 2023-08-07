#' Plot the forest plot for a given meta-analysis
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
#' \code{\link[hMean]{kTRMu}} as the p-value function, and \code{"Pearson"} uses
#' \code{\link[hMean]{pPearsonMu}} as the p-value function.
#' @param heterogeneity The heterogeneity model used to adjust the standard
#' errors \code{se}. Valid options are any combination of \code{"none"},
#' \code{"additive"}, or \code{"multiplicative"}. See also
#' \code{\link[hMean]{hMeanChiSqMu}}.
#' @param diamond_height The height of the polygons. Must be a numeric vector of
#' length 1. Defaults to 0.5.
#' @param studyNames Either \code{NULL} (the default) or a character vector of
#' the same length as \code{thetahat} indicating the names of the individual
#' studies.
#' @param xlim Either \code{NULL} (default) or a numeric vector of length 2,
#' indicating the limits of the x-axis.
#' @param v_space A numeric vector of length 1 indicating the vertical space
#' between the different rows of the plot. Default is 1.5.
#' @param show_studies Must be either \code{TRUE} (default) or \code{FALSE}.
#' If \code{TRUE}, the confidence intervals for the individual studies will
#' be shown in the forest plot. Otherwise, the plot will only show the diamonds
#' for the meta-analyses.
#' @param scale_diamonds Must be either \code{TRUE} or \code{FALSE} with
#' \code{FALSE} being the default. If \code{TRUE}, the maximum height of the
#' diamond is always \code{diamond_height}, even if the p-value function never
#' reaches 1. If \code{FALSE}, the the maximum height of the diamond is
#' \code{diamond_heigth} * max(p-value function).
#'
#' @return An object of class \code{ggplot}. The object contains everything
#' necessary to plot the forest plot.
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
#' ForestPlot(
#'     thetahat = thetahat,
#'     se = se,
#'     distr = distr,
#'     pValueFUN = pValueFUN,
#'     heterogeneity = heterogeneity
#' )
#'
#' @export
ForestPlot <- function(
    thetahat,
    se,
    level = 0.95,
    distr = c("chisq", "f"),
    pValueFUN = c("hMean", "k-Trials", "Pearson", "Edgington", "Fisher"),
    heterogeneity = c("none", "additive", "multiplicative"),
    diamond_height = 0.5,
    v_space = 1.5,
    studyNames = NULL,
    xlim = NULL,
    show_studies = TRUE,
    scale_diamonds = FALSE
) {

    # get the p-value function(s)
    pValueFUN <- match.arg(pValueFUN, several.ok = TRUE)

    # get the desired heterogeneity
    heterogeneity <- match.arg(heterogeneity, several.ok = TRUE)

    # get the distribution
    distr <- match.arg(distr, several.ok = TRUE)

    # make a table with the study intervals
    studyNames  <- if (is.null(studyNames)) paste("Study", seq_along(thetahat))
    se_term <- stats::qnorm(level) * se

    # Make a data frame for the single studies
    studyCIs <- data.frame(
        lower = thetahat - se_term,
        upper = thetahat + se_term,
        estimate = thetahat,
        name = studyNames,
        plottype = 0L,
        color = 0L,
        stringsAsFactors = FALSE,
        row.names = NULL
    )

    # Get the dataframe for the polygons of the new methods
    new_method_cis <- get_CI_new_methods(
        thetahat = thetahat,
        se = se,
        heterogeneity = heterogeneity,
        pValueFUN = pValueFUN,
        distr = distr,
        level = level,
        diamond_height = diamond_height,
        scale_diamonds = scale_diamonds
    )

    # Get the dataframe for the polygons of the old methods
    old_methods_cis <- get_CI_old_methods(
        thetahat = thetahat,
        se = se,
        level = level,
        diamond_height = diamond_height
    )

    # Assemble p_0
    p_0 <- rbind(old_methods_cis$p_0, new_method_cis$p_0)

    # Assemble the dataset
    polygons <- rbind(old_methods_cis$CIs, new_method_cis$CIs)

    # Some cosmetics
    ## If any of the CIs does not exist, replace the row
    na_cis <- anyNA(polygons$x)
    if (na_cis) {
        # Which methods have an NA CI
        na_methods <- with(polygons, unique(name[is.na(x)]))
        # Mark it with a new variable
        polygons$ci_exists <- ifelse(polygons$name %in% na_methods, FALSE, TRUE)
        # Set y coordinate to 0
        idx <- with(polygons, name %in% na_methods)
        polygons$y[idx] <- 0
        polygons$x[idx] <- 0
    } else {
        polygons$ci_exists <- TRUE
    }
    ## Assign correct y-coordinate
    methods <- unique(polygons$name)
    n_methods <- length(methods)
    spacing <- seq(
        (nrow(studyCIs) + n_methods + 1L) * v_space,
        v_space,
        by = -v_space
    )
    studyCIs$y <- spacing[seq_len(nrow(studyCIs))]
    method_spacing <- spacing[(nrow(studyCIs) + 2L):length(spacing)]
    for (i in seq_along(methods)) {
        idx <- with(polygons, name == methods[i])
        polygons$y[idx] <- polygons$y[idx] + method_spacing[i]
    }
    ## Manage colors
    n_colors <- length(unique(polygons$color)) - 1L
    if (n_colors > 0) {
        colors <- scales::hue_pal()(n_colors)
        col_idx <- with(polygons, unique(color[ci_exists & color != 0]))
        colors <- c("gray20", colors[col_idx])
    }
    polygons$color <- factor(polygons$color)

    # Make the plot
    p <- ggplot2::ggplot()
    if (show_studies) {
        p <- p +
            ggplot2::geom_errorbarh(
                data = studyCIs,
                ggplot2::aes(
                    y = y,
                    xmin = lower,
                    xmax = upper,
                    height = diamond_height
                ),
                show.legend = FALSE
            ) +
            ggplot2::geom_point(
                data = studyCIs,
                ggplot2::aes(x = estimate, y = y)
            )
    }
    p <- p +
        ggplot2::geom_polygon(
            data = polygons[polygons$ci_exists == TRUE, ],
            ggplot2::aes(
                x = x,
                y = y,
                group = paste0(name, ".", id),
                fill = color
            ),
            show.legend = FALSE
        ) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title.y = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_line(linetype = "solid"),
            axis.ticks.x = ggplot2::element_line(linetype = "solid"),
            panel.border = ggplot2::element_blank(),
            axis.ticks = ggplot2::element_blank(),
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major.y = ggplot2::element_blank()
        ) +
        ggplot2::scale_y_continuous(
            breaks = spacing,
            labels = c(studyCIs$name, "", unique(polygons$name))
        ) +
        ggplot2::labs(
            x = bquote(mu)
        ) +
        ggplot2::scale_fill_discrete(type = colors)

    if (na_cis) {
        na_rows <- polygons[polygons$ci_exists == FALSE, ]
        for (i in seq_len(nrow(na_rows))) {
            p <- p + ggplot2::annotate(
                geom = "text",
                x = na_rows$x[i],
                y = na_rows$y[i],
                label = "CI does not exist"
            )
        }
    }

    if (!is.null(xlim)) {
        p <- p + ggplot2::xlim(xlim)

    }

    list(
        plot = p,
        p_0 = p_0
    )
}


get_CI_new_methods <- function(
    thetahat,
    se,
    heterogeneity,
    pValueFUN,
    distr,
    level,
    diamond_height,
    scale_diamonds
) {

    # Set up the grid to loop over
    grid <- make_grid(
        pValueFUN = pValueFUN,
        heterogeneity = heterogeneity,
        distr = distr
    )

    # Calculate phi and tau2 (This is the same for all scenarios)
    phi <- hMean::estimatePhi(thetahat, se)
    tau2 <- hMean::estimateTau2(thetahat, se)

    out <- lapply(seq_len(nrow(grid)), function(r) {

        # determine whether to pass phi or tau2
        het <- grid$heterogeneity[r]
        phi_pass <- if (het == "multiplicative") phi else NULL
        tau2_pass <- if (het == "additive") tau2 else NULL
        # determine distr
        distr <- grid$distr[r]
        # get the p-Value function
        pValueFUN <- get(grid$fun_name[r])
        # determine the arguments for the pvalue function
        pValueFUN_args <- list(
            heterogeneity = het,
            phi = phi_pass,
            tau2 = tau2_pass,
            check_inputs = FALSE
        )
        if (!is.na(distr)) {
            pValueFUN_args <- append(pValueFUN_args, list(distr = distr))
        }
        # Calculate the CIs and the minima
        res <- hMeanChiSqCI(
            thetahat = thetahat,
            se = se,
            level = level,
            alternative = "none",
            pValueFUN = pValueFUN,
            pValueFUN_args = pValueFUN_args
        )

        # Calculate polygons
        ci_exists <- if (all(is.na(res$CI))) FALSE else TRUE
        if (ci_exists) {
            # evaluate pValueFUN at thetahat to get the maxima for the diamond
            f_thetahat <- do.call(
                "pValueFUN",
                append(
                    list(thetahat = thetahat, se = se, mu = thetahat),
                    pValueFUN_args
                )
            )
            # Extract the CIs and minima
            CI  <- res$CI
            gamma  <- res$gamma
            # Calculate the polygons for the diamond
            polygons <- calculate_polygon_2(
                CIs = CI,
                thetahat = thetahat,
                f_thetahat = f_thetahat,
                gammas = gamma,
                diamond_height = diamond_height,
                level = level,
                scale_diamonds = scale_diamonds
            )
            polygons$name <- grid$name[r]
            polygons$color <- r

        } else {
            polygons <- data.frame(
                x = NA_real_,
                y = NA_real_,
                id = NA_real_,
                name = grid$name[r],
                color = r
            )
        }

        # Calculate the p-value at mu = 0
        p_0_args <- list(thetahat = thetahat, se = se, mu = 0)
        p_0 <- do.call("pValueFUN", append(p_0_args, pValueFUN_args))
        p_0 <- data.frame(
            name = grid$name[r],
            p_0 = p_0
        )

        # Return
        list(
            CIs = polygons,
            p_0 = p_0
        )
    })

    # Reorganize list
    CIs <- do.call("rbind", lapply(out, "[[", i = 1L))
    p_0 <- do.call("rbind", lapply(out, "[[", i = 2L))

    # Return
    list(CIs = CIs, p_0 = p_0)
}


get_CI_old_methods <- function(
    thetahat,
    se,
    level,
    diamond_height
) {

    # Get the object
    get_obj_reml <- function(thetahat, se, level) {
        meta::metagen(
            TE = thetahat, seTE = se, sm = "MD",
            level = level, method.tau = "REML"
        )
    }
    get_obj_hk <- function(thetahat, se, level) {
        meta::metagen(
            TE = thetahat, seTE = se, sm = "MD",
            level = level, method.tau = "REML", hakn = TRUE
        )
    }
    get_obj_hc <- function(thetahat, se, level) {
        metafor::hc(
            object = metafor::rma(yi = thetahat, sei = se, level = level)
        )
    }
    get_obj <- function(method, thetahat, se, level) {
        switch(
            method,
            "Random Effects, REML" = get_obj_reml(
                thetahat = thetahat,
                se = se,
                level = level
            ),
            "Hartung & Knapp" = get_obj_hk(
                thetahat = thetahat,
                se = se,
                level = level
            ),
            "Henmi & Copas" = get_obj_hc(
                thetahat = thetahat,
                se = se,
                level = level
            )
        )
    }
    # Get the CI
    get_ci_reml <- function(reml) {
        m <- matrix(c(reml$lower.random, reml$upper.random), ncol = 2L)
        colnames(m) <- c("lower", "upper")
        m
    }
    get_ci_hk <- function(hk) {
        m <- matrix(c(hk$lower.random, hk$upper.random), ncol = 2L)
        colnames(m) <- c("lower", "upper")
        m
    }
    get_ci_hc <- function(hc) {
        m <- matrix(c(hc$ci.lb, hc$ci.ub), ncol = 2L)
        colnames(m) <- c("lower", "upper")
        m
    }
    get_ci <- function(method, obj) {
        switch(
            method,
            "Random Effects, REML" = get_ci_reml(reml = obj),
            "Hartung & Knapp" = get_ci_hk(hk = obj),
            "Henmi & Copas" = get_ci_hc(hc = obj)
        )
    }

    # Get the p-value for the null-effect
    get_pval_reml <- get_pval_hk <- function(obj) {
        obj$pval.random
    }
    get_pval_hc <- function(obj, level) {
        ci <- get_ci_hc(obj)
        ReplicationSuccess::ci2p(
            lower = ci[, "lower"],
            upper = ci[, "upper"],
            conf.level = level,
            ratio = FALSE,
            alternative = "two.sided"
        )
    }
    get_pval <- function(method, obj, level) {
        switch(
            method,
            "Random Effects, REML" = get_pval_reml(obj = obj),
            "Hartung & Knapp" = get_pval_hk(obj = obj),
            "Henmi & Copas" = get_pval_hc(obj = obj, level = level)
        )
    }

    # Make a table with the classic methods
    other_methods <- c(
        "Random Effects, REML", "Hartung & Knapp", "Henmi & Copas"
    )

    # Calculate all the measures we want
    ## CIs
    l <- length(other_methods) # The number of methods
    n_pts <- 4L                # The number of points for each CI
    l_tot <- l * n_pts         # The total number of rows
    # Data frame columns for list element CIs
    cis_x <- numeric(l_tot)    # The x coordinates
    cis_y <- numeric(l_tot)    # The y coordinates
    cis_id <- rep(1, l_tot)    # The polygon id
    cis_name <- character(l_tot) # The method name
    cis_color <- rep(0, l_tot) # The polygon color
    # Data frame columns for p_0
    p_0_name <- character(l) # The name column
    p_0_p_0 <- numeric(l) # The p_0 column

    # Calculate stuff for all methodsk
    cis_cnt <- 1L:n_pts
    p_0_cnt <- 1L
    for (meth in other_methods) {
        # Fit the model
        obj <- get_obj(
            method = meth,
            thetahat = thetahat,
            se = se,
            level = level
        )
        # Get the CI and estimate
        ci <- get_ci(method = meth, obj = obj)
        est <- mean(ci)
        # Calculate the p-value from the CI
        p_0_p_0[p_0_cnt] <- get_pval(method = meth, obj = obj, level = level)
        p_0_name[p_0_cnt] <- meth
        # Convert CI to polygon and store in vectors above
        cis_x[cis_cnt] <- c(ci[, "lower"], est, ci[, "upper"], est)
        cis_y[cis_cnt] <- c(0, -diamond_height / 2, 0, diamond_height / 2)
        cis_name[cis_cnt] <- meth
        # Get indices for next iteration
        cis_cnt <- cis_cnt + n_pts
        p_0_cnt <- p_0_cnt + 1L
    }

    list(
        CIs = data.frame(
            x = cis_x,
            y = cis_y,
            id = cis_id,
            name = cis_name,
            color = cis_color,
            stringsAsFactors = FALSE,
            row.names = NULL
        ),
        p_0 = data.frame(
            name = p_0_name,
            p_0 = p_0_p_0,
            stringsAsFactors = FALSE,
            row.names = NULL
        )
    )

    # do.call(
    #     "rbind",
    #     lapply(
    #         other_methods,
    #         function(x) {
    #             obj <- get_obj(
    #                 method = x,
    #                 thetahat = thetahat,
    #                 se = se,
    #                 level = level
    #             )
    #             ci <- get_ci(method = x, obj = obj)
    #             est <- mean(ci)
    #             data.frame(
    #                 x = c(ci[, "lower"], est, ci[, "upper"], est),
    #                 y = c(0, -diamond_height / 2, 0, diamond_height / 2),
    #                 id = 1,
    #                 name = x,
    #                 color = 0
    #             )
    #         }
    #     )
    # )
}

# # Assemble a data frame containing the points for a polygon
calculate_polygon  <- function(CIs, thetahat, gammas = NULL, height, level) {
    # Filter out minima lower than alpha
    if (!is.null(gammas))
        gammas <- gammas[gammas[, 2L] >= (1 - level), , drop = FALSE]
    ## get the x-values
    xVals  <- c(CIs, thetahat)
    if (!is.null(gammas))
        xVals <- c(xVals, gammas[, 1L])
    xVals <- unname(xVals)
    ## precompute some things
    lci <- length(CIs)
    lth <- length(thetahat)
    lx <- length(xVals)
    ## compute the y values
    is_lower <- seq_along(CIs[, 1L])
    is_upper <- is_lower + lci / 2L
    is_theta <- (lci + 1L):(lci + lth)
    is_minimum <- lx - lci - lth
    is_minimum <- if (is_minimum == 0L) {
        integer(0L)
    } else {
        (lx - is_minimum + 1L):lx
    }
    yVals <- double(lx)
    yVals[is_lower | is_upper] <- 0
    yVals[is_theta] <- -height / 2
    if (!is.null(gamma))
        yVals[is_minimum] <- -((gammas[, 2L] - (1 - level)) / level) *
            (height / 2)
    ## reorder vectors according to the x coordinates
    o <- order(xVals)
    xVals <- xVals[o]
    yVals <- yVals[o]
    ## assign polygon IDs based on lower and upper bounds
    is_lower <- which(o %in% is_lower)
    polyID <- integer(lx)
    current_polygon <- 0L
    for (i in seq_along(polyID)) {
        if (i %in% is_lower)
            current_polygon  <- current_polygon + 1L
        polyID[i]  <- current_polygon
    }
    ## repeat parts of the vectors such that we get
    ## the lower and upper points
    poly <- lapply(unique(polyID), function(z) {
        currentX <- xVals[polyID == z]
        currentY  <- yVals[polyID == z]
        l <- length(currentX)
        x <- c(currentX, rev(currentX[-c(1L, l)]))
        y <- c(currentY, rev(currentY[-c(1L, l)] * -1))
        out <- matrix(c(rep(z, length(x)), x, y), ncol = 3L)
        colnames(out) <- c("id", "x", "y")
        out
    })
    as.data.frame(do.call("rbind", poly))
}

# Make the polygon df for classic methods

# Make a df that contains the data for the diamonds
calculate_polygon_2 <- function(
    CIs,
    thetahat,
    f_thetahat,
    gammas,
    diamond_height,
    level,
    scale_diamonds
) {

    # If gamma == NA, there is either one or no CI
    no_ci <- if (all(is.na(CIs))) TRUE else FALSE
    no_gamma <- if (all(is.na(gammas))) TRUE else FALSE

    if (no_ci) {
        return(NULL)
    } else {
        # remove f_thetahat where f_thetahat > alpha
        # because those below alpha are not in CI
        keep <- f_thetahat > 1 - level
        f_thetahat <- f_thetahat[keep]
        thetahat <- thetahat[keep]
        # get lengths
        l_thetahat <- length(f_thetahat)
        l_ci <- nrow(CIs)
        if (!no_gamma) {
            # remove gammas where gammas > alpha
            # because those below alpha are not in CI
            gammas <- gammas[gammas[, 2L] > 1 - level, , drop = FALSE]
            l_gamma <- nrow(gammas)
            # set up x coordinates
            x <- c(CIs[, 1L], gammas[, 1L], thetahat, CIs[, 2L])
            # Get the y-coordinates
            y <- vector("numeric", length(x))
            # type of point
            # 0 = lower ci
            # 1 = minimum
            # 2 = maximum
            # 3 = upper ci
            type <- c(
                rep(0L, l_ci),
                rep(1L, l_gamma),
                rep(2L, l_thetahat),
                rep(3L, l_ci)
            )
            y[type == 1L] <- gammas[, 2L]
        } else {
            # set up x coordinates
            x <- c(CIs[, 1L], thetahat, CIs[, 2L])
            # Get the y-coordinates
            y <- vector("numeric", length(x))
            type <- c(
                rep(0L, l_ci),
                rep(2L, l_thetahat),
                rep(3L, l_ci)
            )
        }
        y[type == 0L | type == 3L] <- 0 # 0 if lower/upper
        y[type == 2L] <- f_thetahat
        # rescale the diamonds if needed (such that the max height is always 1)
        if (scale_diamonds) {
            scale_factor <- 1 / max(y)
            y <- y * scale_factor
        }
        # Order the x & y coordinates
        o <- order(x, decreasing = FALSE)
        x <- x[o]
        y <- y[o]
        type <- type[o]
        # Assign polygon id
        is_upper <- which(type == 3L)
        is_lower <- which(type == 0L)
        id <- vector("numeric", length(x))
        counter <- 1L
        for (i in seq_along(is_upper)) {
            id[is_lower[i]:is_upper[i]] <- counter
            counter <- counter + 1L
        }
        # Arrange this such that it gives back a data.frame
        # that can be used with ggplot
        as.data.frame(
            do.call(
                "rbind",
                lapply(
                    unique(id),
                    function(z) {
                        idx <- id == z
                        x_sub <- x[idx]
                        y_sub <- y[idx] * diamond_height / 2
                        l <- length(x_sub)
                        mat <- matrix(
                            NA_real_,
                            ncol = 3,
                            nrow = 2L * l - 2L
                        )
                        colnames(mat) <- c("x", "y", "id")
                        mat[, 1L] <- c(x_sub, rev(x_sub[-c(1L, l)]))
                        mat[, 2L] <- c(-y_sub, rev(y_sub[-c(1L, l)]))
                        mat[, 3L] <- z
                        mat
                    }
                )
            )
        )
    }
}
