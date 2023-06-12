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
#' @param height The height of the polygons. Must be a numeric vector of length
#' 1. Defaults to 0.5.
#' @param studyNames Either \code{NULL} (the default) or a character vector of 
#' the same length as \code{thetahat} indicating the names of the individual
#' studies.
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
#' pValueFUN <- c("hMean", "k-Trials", "Pearson")
#' distr <- c("chisq")
#' heterogeneity <- "none"
#' ForestPlot(
#'     thetahat = thetahat,
#'     se = se,
#'     distr = distr,
#'     pValueFUN = c("hMean", "k-Trials", "Pearson"),
#'     heterogeneity = NULL
#' )
#' 
#' @export
ForestPlot <- function(
    thetahat,
    se,
    level = 0.95,
    distr = c("chisq", "f"),
    pValueFUN = c("hMean", "k-Trials", "Pearson"),
    heterogeneity = c("none", "additive", "multiplicative"),
    height = 0.5,
    studyNames = NULL
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
    l_theta <- length(thetahat)

    studyCIs <- data.frame(
        lower = thetahat - se_term,
        upper = thetahat + se_term,
        estimate = thetahat,
        name = studyNames,
        plottype = 0L,
        color = 0L,
        #polygon = I(rep(list(NA), l_theta)),
        stringsAsFactors = FALSE,
        row.names = NULL
    )

    # Set up the grid to loop over
    grid <- make_grid(
        pValueFUN = pValueFUN,
        heterogeneity = heterogeneity,
        distr = distr
    )
    grid$pretty_name <- vapply(
        grid$fun_name,
        function(x) {
            switch(
                x,
                "hMeanChiSqMu" = "hMean",
                "kTRMu" = "k-Trials",
                "pPearsonMu" = "Pearson"
            )
        },
        character(1L)
    )
    # hMean <- "hMean" %in% pValueFUN
    # ktrmu <- "k-Trials" %in% pValueFUN
    # if (hMean)
    #     grid_hmean <- expand.grid(
    #         fun_name = "hMean",
    #         heterogeneity = heterogeneity,
    #         distr = distr,
    #         stringsAsFactors = FALSE
    #     )
    # if (ktrmu)
    #     grid_ktrmu <- expand.grid(
    #         fun_name = "k-Trials",
    #         heterogeneity = heterogeneity,
    #         distr = NA_character_,
    #         stringsAsFactors = FALSE
    #     )
    # if (hMean && ktrmu) {
    #     grid <- rbind(grid_hmean, grid_ktrmu)
    # } else if (hMean && !ktrmu) {
    #     grid <- grid_hmean
    # } else {
    #     grid <- grid_ktrmu
    # }
    grid$name <- with(
        grid,
        make_names(
            FUN = pretty_name,
            heterogeneity = heterogeneity,
            distr = distr
        )
    )

    # Calculate phi and tau2
    phi <- hMean::estimatePhi(thetahat, se)
    tau2 <- hMean::estimateTau2(thetahat, se)

    # Calculate the CIs for hMean and k-Trial
    hMeanCIs <- lapply(seq_len(nrow(grid)), function(x) {
        # determine whether to pass phi or tau2
        het <- grid$heterogeneity[x]
        phi_pass <- if (het == "multiplicative") phi else NULL
        tau2_pass <- if (het == "additive") tau2 else NULL
        # get the p-Value function
        pValueFUN <- switch(
            grid$pretty_name[x],
            "hMean" = hMean::hMeanChiSqMu,
            "k-Trials" = hMean::kTRMu,
            "Pearson" = hMean::pPearsonMu
        )
        # Put arguments in a list
        arglist <- list(
            thetahat = thetahat,
            se = se,
            level = level,
            pValueFUN = pValueFUN,
            pValueFUN_args = list(
                heterogeneity = het,
                phi = phi_pass,
                tau2 = tau2_pass
            )
        )
        if (grid$fun_name[x] == "hMean") {
            arglist$pValueFUN_args <- append(
                arglist$pValueFUN_args,
                list(
                    alternative = "none",
                    distr = grid$distr[x]
                )
            )
        }
        # Call hMeanChiSqCI on the arguments
        res <- do.call(hMean::hMeanChiSqCI, arglist)
        # extract the CIs
        CI  <- res$CI
        # Compute the polygons
        gamma  <- res$gamma
        polygons <- calculate_polygon(
            CIs = CI,
            thetahat = thetahat,
            gammas = gamma,
            height = height,
            level = level
        )

        ## Assemble output
        list(
            CI = data.frame(
                lower = rep(CI[, 1L], l_theta),
                upper = rep(CI[, 2L], l_theta),
                estimate = rep(thetahat, each = nrow(res$CI)),
                name = grid$name[x],
                plottype = 1L,
                color = x,
                stringsAsFactors = FALSE,
                row.names = NULL
            ),
            gamma = polygons
        )
    })
    hmeanCIs <- do.call(`rbind`, lapply(hMeanCIs, `[[`, i = 1L))
    hmeanPoly  <- stats::setNames(lapply(hMeanCIs, `[[`, i = 2L), grid$name)
    hmeanPoly <- do.call(
        `rbind`,
        lapply(seq_along(hmeanPoly), function(z) {
            hmeanPoly[[z]]$name <- names(hmeanPoly)[z]
            hmeanPoly[[z]]
        })
    )

    # Calculate CIs for the classic methods
    ## REML
    reml <- meta::metagen(
        TE = thetahat, seTE = se, sm = "MD",
        level = level, method.tau = "REML"
    )
    ## Hartung & Knapp
    hk <- meta::metagen(
        TE = thetahat, seTE = se, sm = "MD",
        level = level, method.tau = "REML", hakn = TRUE
    )
    ## Henmi & Copas
    hc <- metafor::hc(
        object = metafor::rma(yi = thetahat, sei = se, level = level)
    )

    # Make a table with the classic methods
    other_methods <- c(
        "Random Effects, REML", "Hartung & Knapp", "Henmi & Copas"
    )
    other_CIs <- data.frame(
        lower = c(reml$lower.random, hk$lower.random, hc$ci.lb),
        upper = c(reml$upper.random, hk$upper.random, hc$ci.ub),
        estimate = c(
            mean(c(reml$lower.random, reml$upper.random)),
            mean(c(hk$lower.random, hk$upper.random)),
            mean(c(hc$ci.lb, hc$ci.ub))
        ),
        name = other_methods,
        plottype = 1L,
        color = 0L,
        stringsAsFactors = FALSE,
        row.names = NULL
    )
    # Compute the polygons for these methods as well
    otherPoly <- stats::setNames(
        apply(
            as.matrix(other_CIs[c("lower", "upper", "estimate")]),
            1L,
            function(z) {
                dim(z) <- c(1L, 3L)
                out <- calculate_polygon(
                    CIs = z[, 1:2, drop = FALSE],
                    thetahat = z[, 3L],
                    height = height,
                    level = level
                )
                out
            }
        ),
        other_methods
    )
    otherPoly <- do.call(
        `rbind`,
        lapply(seq_along(otherPoly), function(z) {
            otherPoly[[z]]$name <- names(otherPoly)[z]
            otherPoly[[z]]
        })
    )

    # Compose data sets for plotting
    dataCIs <- rbind(
        studyCIs,
        other_CIs,
        hmeanCIs
    )

    dataPoly <- rbind(
        otherPoly,
        hmeanPoly
    )

    # Construct the plot
    ## Calculate positioning of for each of the intervals
    unique_intervals <- unique(dataCIs$name)
    y  <- data.frame(
        name = unique_intervals,
        y = rev(seq_along(unique_intervals)),
        stringsAsFactors = FALSE
    )
    y$y[seq_along(thetahat)] <- y$y[seq_along(thetahat)] + 1L
    dataCIs <- merge(
        dataCIs, y,
        by = "name",
        sort = FALSE,
        all.x = TRUE,
        all.y = FALSE
    )
    dataPoly <- merge(
        dataPoly, y,
        by = "name",
        sort = FALSE,
        all.x = TRUE,
        all.y = FALSE
    )
    dataPoly$y <- dataPoly$y.x + dataPoly$y.y
    dataPoly <- dataPoly[c("name", "id", "x", "y")]

    # Construct id column for each polygon
    dataPoly$ID <- with(dataPoly, paste(name, id))

    # Construct y-axis labels
    meth <- unique(dataCIs$name)
    ylabs <- rev(c(
        meth[seq_along(thetahat)],
        "",
        meth[(l_theta + 1L):length(meth)]
    ))

    # Add color scales to polygons
    dataPoly <- merge(
        dataPoly,
        dataCIs[c("name", "color")],
        by = "name",
        all.x = TRUE,
        all.y = FALSE,
        sort = FALSE
    )
    dataPoly$color <- as.factor(dataPoly$color)

    # Which CIs to draw as polygons and which as bars
    draw_poly <- unique(dataCIs$name[dataCIs$plottype == 1L])
    draw_bar <- unique(dataCIs$name[dataCIs$plottype == 0L])

    # Do the plot
    ggplot2::ggplot() +
    ggplot2::geom_errorbarh(
        data = dataCIs[dataCIs$name %in% draw_bar, ],
        ggplot2::aes(y = y, xmin = lower, xmax = upper, height = height),
        show.legend = FALSE
    ) +
    ggplot2::geom_point(
        data = dataCIs[dataCIs$name %in% draw_bar, ],
        ggplot2::aes(x = estimate, y = y)
    ) +
    ggplot2::geom_polygon(
        data = dataPoly[dataPoly$name %in% draw_poly, ],
        ggplot2::aes(x = x, y = y, group = ID, fill = color),
        show.legend = FALSE
    ) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::scale_y_continuous(
        breaks = seq(1L, max(dataCIs$y), 1L),
        labels = ylabs
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
        x = bquote(mu)
    ) +
    ggplot2::scale_fill_discrete(
        type = c(
            "gray20",
            scales::hue_pal()(length(unique(dataPoly$color)) - 1L)
        )
    ) +
    ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        panel.border = ggplot2::element_blank(),
        axis.ticks = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank()
    )
}

make_names <- function(FUN, heterogeneity = NULL, distr = NULL) {
    # handle heterogeneity
    if (!is.null(heterogeneity))
        heterogeneity <- vapply(heterogeneity, function(x) {
            switch(
                x,
                "none" = " none",
                "additive" = " add.",
                "multiplicative" = " mult."
            )
        }, character(1L))
    else
        heterogeneity <- rep("", length(FUN))
    # handle distr
    if (!is.null(distr))
        distr <- ifelse(is.na(distr), "", paste0(" (", distr, ")"))
    else
        distr <- rep("", length(FUN))
    # Make names
    paste0(FUN, heterogeneity, distr)
}

# Assemble a data frame containing the points for a polygon
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
    as.data.frame(do.call(`rbind`, poly))
}
