#' @title Visualizations of \code{confMeta} objects
#'
#' @description Plots one or more \code{confMeta} objects. This function can
#'     create two types of plots: the \emph{p}-value function plot (also known
#'     as drapery plot) and the forest plot. This allows for direct visual
#'     comparison of different \emph{p}-value functions.
#'
#' Optionally, a Bayesian meta-analysis object created with the
#' \code{bayesmeta::bayesmeta()} function can be supplied via the \code{bayesmeta} argument.
#' When provided, its posterior summary is displayed as an additional diamond
#' at the bottom of the forest plot for comparison with the frequentist methods
#' 
#' \strong{Important Note:} If supplying a Bayesian meta-analysis object 
#' via the \code{bayesmeta} argument, this function explicitly extracts the 
#' 95% credible intervals. If the Bayesian model was fit using a different 
#' credible level, the function will crash.
#' 
#' @param ... One or more objects of class \code{confMeta}.
#' @param type A character vector of length 1 or 2. Indicates what type of
#'     plot should be returned. Accepted values are \code{"p"}, \code{"forest"},
#'     or both. Defaults to \code{c("p", "forest")}.
#' @param diamond_height Numeric scalar. Indicates the maximal
#'     possible height of the diamonds in the forest plot. Defaults to 0.5.
#'     This argument is only relevant if \code{type} contains \code{"forest"} and will
#'     be ignored otherwise.
#' @param v_space Numeric scalar. Indicates the vertical space
#'     between two diamonds in the forest plot. Defaults to 1.5.
#'     This argument is only relevant if \code{type} contains \code{"forest"} and will
#'     be ignored otherwise.
#' @param scale_diamonds Logical. If \code{TRUE} (default), the diamond is rescaled to the
#'      interval \[0, 1\] in cases where the maximum of the p-value
#'     function is not equal to 1. This argument is only relevant if \code{type}
#'     contains \code{"forest"} and will be ignored otherwise.
#' @param show_studies Logical. If \code{TRUE} (default), the forest plot shows the
#'     confidence intervals for the individual effect estimates. Otherwise,
#'     the intervals are suppressed. This argument is only relevant if \code{type}
#'     contains \code{"forest"} and will be ignored otherwise.
#' @param drapery Logical. If \code{TRUE} (default), individual
#'     study effects are represented as drapery plots in the \emph{p}-value function plot. 
#'     If \code{FALSE} the studies are represented by a simple vertical line at their effect estimates.
#'      This argument is only relevant if \code{type}
#'     contains \code{"p"} and will be ignored otherwise.
#' @param reference_methods_p Character vector controlling which reference
#'     meta-analysis methods are shown as baseline curves in the
#'     \emph{p}-value function plot. Valid options are \code{"fe"},
#'     \code{"re"}, or \code{c("fe", "re")} for both. Defaults to \code{"re"}.
#' @param reference_methods_forest Character vector of length 1 to 4.
#'     Specifies which reference meta-analysis methods should be shown as
#'     diamonds in the forest plot. Valid options are any subset of
#'     \code{c("fe", "re", "hk", "hc")}, which correspond to:
#'     \itemize{
#'       \item \code{"fe"}: fixed-effect meta-analysis
#'       \item \code{"re"}: random-effects meta-analysis
#'       \item \code{"hk"}: Hartung-Knapp adjustment
#'       \item \code{"hc"}: Henmi-Copas adjustment
#'     }
#'     Defaults to \code{c("re", "hk", "hc")}.
#' @param ref_labels A named character vector to customize the labels of the reference 
#'     methods in the plot legends and axes (e.g., \code{c("fe" = "Fixed-effect (IV)", 
#'     "re" = "Random-effects (DL)")}). Defaults to \code{NULL}.
#' @param xlim Numeric vector of length 2. Global limits for the x-axis. 
#'   If \code{NULL}, limits are calculated automatically.
#' @param xlim_p Numeric vector of length 2. Specific x-axis limits for the 
#'   \emph{p}-value plot. Overrides \code{xlim}.
#' @param xlim_forest Numeric vector of length 2. Specific x-axis limits for 
#'   the forest plot. Overrides \code{xlim}.
#' @param same_xlim Logical. If \code{TRUE}, forces the \emph{p}-value plot to use the 
#'   same x-axis limits as the forest plot. Defaults to \code{TRUE}.
#' @param xlab Character string. Label for the x-axis. Defaults to \code{NULL} 
#'   (renders as \eqn{\mu}).
#' @param n_breaks Numeric. Approximate number of tick marks to display 
#'   on the x-axis. Defaults to 7.
#' @param bayesmeta An object of class \code{bayesmeta}, typically created using 
#'     \code{bayesmeta::bayesmeta()}. When provided, the posterior median (or mean) and
#'     95% credible interval are displayed as an additional diamond at the bottom 
#'     of the forest plot. Has no effect on the \emph{p}-value plot. Defaults to \code{NULL}.
#' @param n_points Numeric. Number of points used to create the \emph{p}-value plot. 
#'     Higher values (e.g., 10000) yield higher resolution but take longer to render. 
#'     Defaults to 1000.
#' @return An object of class \code{ggplot} containing the specified plot(s).
#'
#'
#'
#'
#' @importFrom patchwork wrap_plots
#'
#' @export
autoplot.confMeta <- function(
    ...,
    type = c("p", "forest"),
    diamond_height = 0.5,
    v_space = 1.5,
    scale_diamonds = TRUE,
    show_studies = TRUE,
    drapery = TRUE,
    # drapery_ma = c("fe"),
    reference_methods_p      = c("re"),
    reference_methods_forest = c("re", "hk", "hc"),
    ref_labels = NULL,
    show_unadjusted = TRUE,
    xlim = NULL,
    xlim_p = NULL,
    xlim_forest = NULL,
    same_xlim = TRUE,
    xlab = NULL,
    n_breaks = 7,
    bayesmeta = NULL,
    n_points = 1000
) {

    # Check validity of reference methods
  check_ref_methods(reference_methods = reference_methods_forest)
  if (!all(reference_methods_p %in% c("fe", "re"))) {
    stop("`reference_methods_p` can only contain \"fe\", \"re\", or both.")
  }
  type <- match.arg(type, several.ok = TRUE)
  reference_methods_forest <- match.arg(
    reference_methods_forest,
    several.ok = TRUE,
    choices = c("fe", "re", "hk", "hc")
  )
  reference_methods_p <- match.arg(
    reference_methods_p,
    several.ok = TRUE,
    choices = c("fe", "re")
  )
  
  
  # Check for ref_labels
  valid_method_keys <- c("fe", "re", "hk", "hc")
  
  if (!is.null(ref_labels)) { 
    if (!is.character(ref_labels) || is.null (names(ref_labels))) {
      stop("'ref_labels' must be a named character vector, e.g. c(fe = 'Fixed-effect (IV)', re = 'Random-effects (REML)')")
    }
    invalid_keys <- setdiff(names(ref_labels), valid_method_keys)
    if (length(invalid_keys) > 0) {
      stop(paste0("Invalid names in `ref_labels`: ", paste(invalid_keys, collapse = ", "),
                  ". Must be one of: ", paste(valid_method_keys, collapse = ", "), "."))
    }
  }
  

    # Check all the confMeta objects
    cms <- list(...)
    check_cms(cms = cms)
    check_unique_names(cms = cms)

    # Check that the levels, estimates, SEs are equal
    # in all confMeta objects
    ests <- lapply(cms, "[[", i = "estimates")
    ses <- lapply(cms, "[[", i = "SEs")
    lvl <- lapply(cms, "[[", i = "conf_level")
    nms <- lapply(cms, "[[", i = "study_names")
    ks <- lapply(cms, "[[", i = "k_studies")
    ok <- c(
        check_equal(ests),
        check_equal(ses),
        check_equal(lvl),
        check_equal(nms),
        check_equal(ks)
    )
    if (any(!ok)) {
        stop(
            paste0(
                "For plotting, all confMeta objects must have the same ",
                "'estimates', 'SEs', 'study_names', 'conf_level' and 'k_studies' elements."
            )
        )
    }

    # Add some more checks for other graphics parameters
    check_TF(x = scale_diamonds)
    check_TF(x = show_studies)
    check_TF(x = drapery)
    check_length_1(x = diamond_height)
    check_class(x = diamond_height, "numeric")
    check_length_1(x = v_space)
    check_class(x = v_space, "numeric")

    # if xlim not specified, builds a range that covers all CIs and add 5% of margin
    # I added the possibility to handle xlim_p and xlim_forest divvertently!
    
    # COLOR HANDLING
    fun_names_shared  <- vapply(cms, "[[", i = "fun_name", character(1L)) #take the fun_name from every confMeta object (Edgington, Edgington_w ...)
    ref_base_p        <- intersect(c("fe", "re"), reference_methods_p)
    ref_forest_only   <- setdiff(
      map_ref_methods(reference_methods_forest, labels = ref_labels),
      map_ref_methods(ref_base_p, labels = ref_labels)
    )
    fac_levels_all    <- c(map_ref_methods(ref_base_p, labels = ref_labels), ref_forest_only, fun_names_shared)
    shared_colors     <- c("#000000", scales::hue_pal()(length(fac_levels_all) - 1))
    names(shared_colors) <- fac_levels_all #attach each method name to the color vector
    
    #This is the complete list of everything that needs a color 
    #we want one color x method -> first is black, then hue_pal() generates evenly sapced color
    
    
    # xlim for p plot:
    if (!is.null(xlim_p)) {
        check_xlim(x = xlim_p)
    } else if (!is.null(xlim)) { #goes to the global xlim if specified
      xlim_p <- xlim
      check_xlim(xlim_p)
    } else {
      candidates_p <- unname(
            do.call(
                "c",
                lapply(cms, function(x) {
                  
                  comp <- x$comparison_cis
                  
                  # For the p value plot, we only care about FE or RE. So we don't consider the other comparison method
                  keep <- rownames(comp) %in% c("Fixed-effect", "Random-effects", "fe", "re")
                  comp_filtered <- comp[keep, , drop = FALSE]
                  
                  with(x, c(individual_cis, joint_cis, comp_filtered))
                    }
                )))
        
        candidates_p <- candidates_p[is.finite(candidates_p)]
        
        ext_perc <- 5
        lower_p <- min(candidates_p)
        upper_p <- max(candidates_p)
        margin_p <- (upper_p - lower_p) * ext_perc / 100
        xlim_p <- c(lower_p - margin_p, upper_p + margin_p)
    }
    
    
    # xlim for f plot:
    if (!is.null(xlim_forest)) {
      
      check_xlim(xlim_forest) 
    } else if (!is.null(xlim)) {
      xlim_forest <- xlim     
      check_xlim(xlim_forest)
      
    } else {
      candidates_f <- unname(do.call("c", lapply(cms, function(x) {
        
        comp <- x$comparison_cis
        
        # create x axis only on the actual ref methods
        allowed_methods <- map_ref_methods (reference_methods_forest)
        
        keep <- rownames(comp)%in%allowed_methods
        comp_filtered <- comp[keep, , drop = FALSE]
        
        
        with(x, c(individual_cis, joint_cis, comp_filtered))
      })))
      
      candidates_f <- candidates_f[is.finite(candidates_f)]
      
      if (length(candidates_f) == 0) {
        xlim_forest <- c(-1, 1)
      } else {
        diff_f <- max(candidates_f) - min(candidates_f)
        margin_f <- if(diff_f == 0) 0.5 else diff_f * 0.05
        xlim_forest <- c(min(candidates_f) - margin_f, max(candidates_f) + margin_f)
      }
    }
    
    
    # if you want to have same x lim, set the x lim_p (which can be different since we excluded some methods) = lim_f
    if (same_xlim) {
      xlim_p = xlim_forest
    }

    plots <- list()
    if ("p" %in% type) {
      plots[["p"]] <- ggPvalueFunction(
        cms=cms, 
        drapery = drapery, 
        xlim = xlim_p,
        xlab = xlab,
        reference_methods = reference_methods_p, 
        n_breaks =  n_breaks,
        n_points = n_points,
        shared_colors = shared_colors,
        ref_labels = ref_labels
        
        
      )
    }
    
    if ("forest" %in% type) {
      plots[["forest"]] <- ForestPlot(
        cms = cms, 
        diamond_height = diamond_height, 
        v_space = v_space, 
        xlim = xlim_forest, 
        show_studies = show_studies, 
        scale_diamonds = scale_diamonds, 
        reference_methods = reference_methods_forest, 
        xlab = xlab, 
        n_breaks =  n_breaks,
        shared_colors = shared_colors,
        ref_labels = ref_labels,
        show_unadjusted = show_unadjusted
        
        
      )
    }
    
    
    # Add Bayesmeta if requested
    if (!is.null(bayesmeta)) {
      if (!inherits(bayesmeta, "bayesmeta")) stop("bayesmeta argument MUST be of class 'bayesmeta'")
           if ("forest" %in% type) {
             plots[["forest"]] <- add_bayes_forest(
               p = plots[["forest"]], 
               bm = bayesmeta, 
               diamond_height = diamond_height,
               v_space = v_space,               
               color = "black", 
               label = "Bayesmeta"
             )
                }
    }
    
    # OUTPUT
    if (length(plots) == 0) {
      return(NULL)
    } else if (length(plots) == 1L) {
      return(plots[[1L]])
    } else {
      return(patchwork::wrap_plots(plots, ncol = 1L)) #put them vertically if there are 2 plots
    }
}
    

# ==============================================================================
# Reexport autoplot function
# ==============================================================================

#' @importFrom ggplot2 autoplot
#' @export
ggplot2::autoplot

# ==============================================================================
# Input checks
# ==============================================================================

check_cms <- function(cms) {
  for (i in seq_along(cms)) {
    msg <- tryCatch(
      validate_confMeta(cms[[i]]),
      error = function(e) conditionMessage(e)
    )
    if (!is.null(msg)) {
      stop(
        paste0("Problem found in confMeta object number ", i, ": ", msg),
        call. = FALSE
      )
    }
  }
  invisible(NULL)
}

check_unique_names <- function(cms) {
    nms <- vapply(cms, "[[", i = "fun_name", character(1L))
    if (any(duplicated(nms))) {
        stop(
            "confMeta objects must have unique `fun_name`.",
            call. = FALSE
        )
    }
    invisible(NULL)
}

check_equal <- function(x) {
  x    <- lapply(x, sort, decreasing = FALSE)
  comp <- x[[1L]]
  all(vapply(x, function(xi) all(xi == comp), logical(1L)))
}

check_TF <- function(x) {
    obj <- deparse1(substitute(x))
    if (!is_TF(x)) {
        stop(
            paste0(
                "Object ", obj, " must be either `TRUE` or `FALSE`."
            ),
            call. = FALSE
        )
    }
    invisible(FALSE)
}

is_TF <- function(x) {
    if (length(x) == 1L && is.logical(x) && !is.na(x)) {
        TRUE
    } else {
        FALSE
    }
}

check_xlim <- function(x) {
    ok <- length(x) == 2L && is.numeric(x) && x[1L] < x[2L]
    if (!ok) {
        stop(
            paste0(
                "Argument `xlim` must be a numeric vector of length 2 where ",
                "the first element is smaller than the second."
            ),
            call. = FALSE
        )
    }
    invisible(NULL)
}


check_ref_methods <- function(reference_methods) {
    # Check for invalid reference_methods
    valid_methods <- c("fe", "re", "hk", "hc")
    is_invalid <- !(reference_methods %in% valid_methods)
    if (any(is_invalid)) {
        msg <- paste0(
            "Detected invalid reference_methods: ",
            paste0(reference_methods[is_invalid], collapse = ", "),
            "."
        )
        stop(msg)
    }
}

# ==============================================================================
# p-value function plot
# ==============================================================================

#' @importFrom stats qnorm
#' @importFrom ggplot2 ggplot aes xlim geom_line geom_hline geom_vline
#' @importFrom ggplot2 geom_point scale_y_continuous sec_axis geom_segment
#' @importFrom ggplot2 theme theme_minimal labs element_text element_blank
#' @importFrom ggplot2 xlim scale_color_manual
#' @importFrom scales hue_pal
#' @importFrom stats setNames
#' @noRd
ggPvalueFunction <- function(
    cms,
    xlim,
    drapery,
    xlab,
    reference_methods,
     n_breaks, 
    n_points,
    shared_colors,
    ref_labels
) {
 
  # Set some constants that are equal for all grid rows
  #use const as a list for containing stuff, muSeq is for creating the grid of mu values,

  const <- list(
    estimates  = cms[[1L]]$estimates,
    SEs        = cms[[1L]]$SEs,
    conf_level = cms[[1L]]$conf_level,
    eps = 0.0025,
    eb_height = 0.025,
    muSeq = seq(xlim[1], xlim[2], length.out = n_points)
  )
  
  

    # Get the function names (for legend)
    fun_names <- vapply(cms, "[[", i = "fun_name", character(1L))
    
    ref_base <- intersect(c("fe", "re"), reference_methods)
    fac_levels <- c(map_ref_methods(ref_base, labels = ref_labels), fun_names)
    
    # Calculate the p-values and CIs by looping along each cm object
    data <- lapply(seq_along(cms), function(x) {
        
        cm <- cms[[x]]
        fun <- cm$p_fun
        fun_name <- cm$fun_name
        alpha <- 1 - const$conf_level
        w <- if (!is.null(cm$w)) cm$w else NULL
        k <- cm$k_studies
        
        
        
        # use the new call for the function 
        pval <- call_pfun(
          fun = fun,
          estimates = cm$estimates,
          SEs = cm$SEs,
          mu = const$muSeq,
          w = w,
          k = k 
        )
        
        
        ################################################################################################################
        #         ######## NOTE FOR WHO CARES 
        # (!!) here we call just those arguments, cause the others are already saved inside p_fun since it is a closure 
        #    thanks to make_p_fun() (e.g. heterogeneity, approx, tau2, neff...)
        #
        #############################################################################################à##################
        
        
        
        # now we have a vector p val of length 10.000 (one for each mu) 
        
        
        pval <- pmin(pmax(pval, 0), 1) #sometimes due to approx error are a tiny bit smaller than 0, fix to not have any error messages
        
        CIs <- cm$joint_cis
        y0 <- cm$p_0[, 2L] #p_val at mu = 0

        # Data frame with the lines of the p-value function: (x,y coords, value
        # at x = 0)
        df1 <- data.frame(
            x = const$muSeq,
            y = pval,
            y0 = y0, # p value at mu = 0
            group = factor(fun_name, levels = fac_levels), #to plot different lines each method
            stringsAsFactors = FALSE,
            row.names = NULL
        )
        # make a second data frame for the display of confidence intervals
        df2 <- data.frame(
            xmin = CIs[, 1L],
            xmax = CIs[, 2L],
            # y = rep(1 - conf_level, nrow(CIs$CI)) + factor * const$eps,
            y = rep(alpha, nrow(CIs)), #horiz segment at confidence level
            group = factor(fun_name, levels = fac_levels),
            stringsAsFactors = FALSE,
            row.names = NULL
        )
        #define the height of the small vertical tiks at the extreme of CI:
        #  |
        #--|---------|--
        #  |
        df2$ymax <- df2$y + const$eb_height 
        df2$ymin <- df2$y - const$eb_height
        list(df1, df2)
    })
    
    #now we have the data object (e.g. data[[1]] = list(df1, df2) is the first ConfMeta object
    # data[[2]] = list(df1, df2) second ConfMeta obj etc...)
    
    # Create the whole data frame for the lines with p-value functions
    # as well as the data frame for the error bars
    lines     <- do.call("rbind", lapply(data, "[[", 1L))
    errorbars <- do.call("rbind", lapply(data, "[[", 2L))

    
    
    ## if you print errorbars:
    #       xmin     xmax           y    group       ymax  ymin
    # -0.09517928 1.032893        0.05  Fisher       0.075 0.025 
    #  0.28089910 1.148117       0.05  Edgington     0.075 0.025 

    # means that the graph for Fisher will be drawed: an horiz segm from mu=-0.095 to mu=1.03, positioned at level y = 0.05, with two vertical ticks from 0.025 to 0.075

    # Calculate the drapery lines if TRUE (by default)
    # dp is just a df with two columns: x, y and the study name (if names are nto provided, just the number of the study)
    if (drapery) {
        dp <- get_drapery_df(
            estimates = const$estimates,
            SEs = const$SEs,
            mu = const$muSeq,
            k = cms[[1L]]$k_studies
            
        )
        
        
        ## also compute the p-value function for the comparison methods
        
        ### Object containing CIs for the reference methods
        #NOTE--> this is the exact reason why we must use the same k for all the cms object we input
        # cause otw the comparison methods plotted are the ones of the first object. 
        
        # IN GENERAL REMEMBER --> the comparison results CANNOT vary between the cms provided to the plotting function 
        #they must be all equal!
        
        ci_obj <- cms[[1L]]$comparison_cis   
        cl <- cms[[1L]]$conf_level
        
        ## What reference method should we calculate thedrapery plot for:
        ## "fe" or "re"
        ref_names <- intersect(
          map_ref_methods(c("fe", "re")),          # all possible default names
          map_ref_methods(reference_methods)       # only those requested
        )
        ref_names <- ref_names[ref_names %in% rownames(ci_obj)]
        
        # Map from default name -> custom label (for legend)
        default_to_label <- setNames(
          map_ref_methods(c("fe", "re", "hk", "hc"), labels = ref_labels),
          map_ref_methods(c("fe", "re", "hk", "hc"))
        )
        

        
        

          # since we only have lower and upper in ci_obj, we extract from it the method we want
          # e.g. fe or re. Then compute the mean of lower and upper as the estimate, and reverse engegneer the se. 
          
          # then use get_drapery_df 
        
        rmadf <- do.call("rbind", lapply(ref_names, function(nm) {
          rmaestimate <- mean(ci_obj[nm, ])
          rmase       <- (ci_obj[nm, 2L] - ci_obj[nm, 1L]) /
            (2 * stats::qnorm((1 + cl) / 2))
          df       <- get_drapery_df(estimates = rmaestimate, SEs = rmase, mu = const$muSeq)
          df$study <- NULL
          df$group <- default_to_label[nm]
          
          df
        }))
        rmadf$group <- factor(rmadf$group, levels = fac_levels)
        
    }
    
    # now dp and rmdaf are two similar df, one for my methods, the other for the baseline method  

    # Define function to convert breaks from primary y-axis to
    # breaks for secondary y-axis
    trans <- function(x) abs(x - 1) * 100 #trasform p value p=0.05 -> 95%.
    
    # Define breaks for the primary y-axis (the one on the right)
    b_points <- c(1 - const$conf_level, pretty(c(lines$y, 1)))
    o <- order(b_points, decreasing = FALSE)
    breaks_y1 <- b_points[o]
    
    # Compute breaks for the secondary y-axis
    # breaks_y2 <- trans(b_points[o])
    # Set transparency
    transparency <- 1

    #add plots for single studies
    p <- ggplot2::ggplot(
      data = lines,
      ggplot2::aes(x = x, y = y, color = group)
    ) + 
      ggplot2::geom_hline(yintercept = 1 - const$conf_level, linetype = "dashed") + # significance treshold
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") #line at mu = 0
    
    
    # DEFINE THE axis scale!

    p <- p +
      ggplot2::scale_x_continuous(
        breaks = scales::pretty_breaks(n = 7),  # 7 number is just arbitrary
        minor_breaks = NULL
      ) +
      ggplot2::coord_cartesian(xlim = xlim, ylim = c(0, 1), expand = FALSE) # this allows us to ZOOM in the plot, not cut it 
    
    #now add the plot for baseline method
    if (!drapery) { #only vertical lines if drapery FALSE
        p <- p + ggplot2::geom_vline(
            xintercept = estimates, #or const$estimates
            linetype = "dashed"
        )
    } else { #otwise add the drapery
        p <- p +
            ggplot2::geom_line(
                data = dp,
                mapping = ggplot2::aes(x = x, y = y, group = study),
                linetype = "dashed",
                color = "lightgrey",
                show.legend = FALSE
            ) +
            ggplot2::geom_line(
                data = rmadf,
                mapping = ggplot2::aes(x = x, y = y, color = group),
                # color = "#00000099",
                show.legend = TRUE
            )


    }
    if (0 > xlim[1L] && 0 < xlim[2L]) {
        # Vertical line at 0
        p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
            ggplot2::geom_point(
                # Points at (x = 0, y = p_0) (it is the H0, put a point there)
                data = lines[!is.na(lines$y0), ],
                ggplot2::aes(x = 0, y = y0, color = group),
                alpha = transparency
            )
    }
    p <- p +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") + #horiz line at y=0
        ggplot2::geom_line(alpha = transparency) +
        ggplot2::scale_y_continuous(
            name = "p-value", #left axis 
            breaks = breaks_y1,
            #limits = c(0, 1), now there should be coord_cartesian
            #expand = c(0, 0),
            sec.axis = ggplot2::sec_axis(
                transform = trans,
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
        # ggplot2::labs(
        #     x = bquote(mu),
        #     color = "Configuration"
        # ) +
        # Set theme
        ggplot2::theme_minimal() +
        ggplot2::theme(
            axis.title.y.right = ggplot2::element_text(angle = 90),
            legend.position = "bottom"
        )

    # Add axis/legend label
    xlab <- if (is.null(xlab)) {
        bquote(mu)
    } else {
        xlab
    }
    p <- p +
        ggplot2::labs(
          x = xlab,
          color = NULL
        ) +
        ggplot2::theme(
          legend.title = ggplot2::element_blank()
        )
    

    # Add SHARED coloring for curves
    p <- p + ggplot2::scale_color_manual(values = shared_colors)
    

    
    
    
    # return
    p
}

# Calculate the drapery lines, i.e. p-value functions of single studies, that are plotted as 
# grey dashed lines in the back
#' @importFrom stats pnorm
#' @noRd
get_drapery_df <- function(estimates, SEs, mu, k = rep(1, length(estimates) ) ) {

      # get lenghts
    l_t <- length(estimates) #n studies
    l_m <- length(mu) # n point on grid (e.g. 1000 or 10000)
    l_tot <- l_t * l_m # total n points to generate
    
    # Initialize vectors
    x_dp <- rep(mu, times = l_t)  #create the grids [mu1,...,mu10000, mu1, ..., mu10000, ..] repetead for each study
    y_dp <- study <- numeric(l_tot)
    
    # Indices to loop over
    idx <- seq_len(l_m) # from 1 to 1000
    nms <- names(estimates)
    
    # Loop over each study
    for (i in seq_along(estimates)) {
      
      # Raw two-sided p-value curve for study i
      p_raw <- 2 * (1 - stats::pnorm(abs(estimates[i] - mu) / SEs[i])) # creates 1000 p values
      
      # Apply best-of-k adjustment (identity when k[i] == 1)
      y_dp[idx] <- 1 - (1 - p_raw)^k[i] # in the y_dp vector, which is of length 1000* (number of studies), put the computed p values
      
      study[idx] <- if (is.null(nms)) rep(i, l_m) else rep(nms[i], l_m)
      idx <- idx + l_m  # move to the next indeces, i.e. at the start idx is from 1 to 1000, then becomes from 1001 to 2000
    }
    #combine studies in a single df
    data.frame( 
        x = x_dp,
        y = y_dp,
        study = study,
        stringsAsFactors = FALSE
    )
}

map_ref_methods <- function(abbrevs, labels = NULL) {
  defaults <- c(
    re = "Random-effects",
    hk = "Hartung & Knapp",
    hc = "Henmi & Copas",
    fe = "Fixed-effect"
  )
  
  # UPDATED FUNCTION: drop the switch and set the default dictionary, then if we have labels provided
  # replace them with the names we have
  if (!is.null(labels)) {
    defaults[names(labels)] <- labels
  }
  unname(defaults[abbrevs])
  
}

# ==============================================================================
# forest plot
# ==============================================================================

#' @importFrom ggplot2 ggplot xlim geom_vline geom_errorbar aes geom_point
#' @importFrom ggplot2 geom_polygon theme_minimal theme element_blank
#' @importFrom ggplot2 element_line scale_y_continuous labs scale_fill_discrete
#' @importFrom ggplot2 annotate
#' @importFrom scales hue_pal
#' @importFrom rlang .data
#' @noRd

ForestPlot <- function(
    cms,
    diamond_height,
    v_space,
    show_studies,
    scale_diamonds,
    reference_methods,
    xlim,
    xlab,
    n_breaks,
    shared_colors,
    ref_labels,
    show_unadjusted 
) {

  
     # # Make a data frame for the single studies
     # Also here we make from just the first object, since all cms object must have the same individual cis
     # cannot provide different k !!
    cm <- cms[[1L]]
    z_crit <- stats::qnorm(1 - (1 - cm$conf_level) / 2)
    
    studyCIs <- data.frame(
        lower = cm$individual_cis[, 1L],
        upper = cm$individual_cis[, 2L],
        lower_unadj       = cm$estimates - z_crit * cm$SEs,  # Wald lower
        upper_unadj       = cm$estimates + z_crit * cm$SEs,  # Wald upper
        estimate          = cm$adjusted_estimates,  # bias-corrected estimate
        estimate_raw      = cm$estimates,            # original value
        name = cm$study_names,
        plottype = 0L,
        color = 0L,
        stringsAsFactors = FALSE,
        row.names = NULL
    )

    # Get the dataframe for the polygons of the new methods
    new_method_cis <- get_CI_new_methods(
        cms = cms,
        diamond_height = diamond_height,
        scale_diamonds = scale_diamonds
    )

    # RETRIEVE (not recompute) the dataframe for the polygons of the old methods
    old_methods_cis <- get_CI_old_methods(
      cm = cm,
      diamond_height = diamond_height
    )



    # Rename the default names to ref_labels!
    # first creates named vector (like Fixed effect --> Fixed effect (IV))
    default_to_custom <- setNames(
      map_ref_methods(c("fe", "re", "hk", "hc"), labels = ref_labels),
      map_ref_methods(c("fe", "re", "hk", "hc"))
    )
    
    #replace
    old_methods_cis$CIs$name <- default_to_custom[old_methods_cis$CIs$name] # ? maybe wrap those i unname? 
    old_methods_cis$p_0$name <- default_to_custom[old_methods_cis$p_0$name]
    
    
    reference_methods <- map_ref_methods(abbrevs = reference_methods, labels = ref_labels)

    keep_cis <- with(old_methods_cis, CIs$name %in% reference_methods)
    keep_p0 <- with(old_methods_cis, p_0$name %in% reference_methods)

    old_methods_cis$CIs <- old_methods_cis$CIs[keep_cis, ]
    old_methods_cis$p_0 <- old_methods_cis$p_0[keep_p0, ]

    ###########################################################################
    # Quick fix, remove the above at some point                               #
    ###########################################################################

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
    polygons$color_hex <- ifelse(
      polygons$name %in% names(shared_colors),
      shared_colors[polygons$name],
      "gray20"
    )
    
    # Make the plot
    p <- ggplot2::ggplot()
    
    
    # DEFINE THE axis scale! smartly 
    
    if (!is.null(xlim)) {
      p <- p +
        ggplot2::scale_x_continuous(
          # limits = xlim,
          breaks = scales::pretty_breaks(n = n_breaks),
          minor_breaks = NULL
        ) + 
        ggplot2::coord_cartesian(xlim = xlim, clip = "on", expand = FALSE) # 
      
      # vertical line at 0
      if (0 > xlim[1L] && 0 < xlim[2L]) {
        p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed")
      }
    } else {
      p <- p +
        ggplot2::scale_x_continuous(
          breaks = scales::pretty_breaks(n = n_breaks),
          minor_breaks = NULL
        )
    }
    
    if (show_studies) {
      if (show_unadjusted && any(cm$k_studies != 1)) {
        
        y_offset <- diamond_height * 0.7   # vertical shift to avoid overlap
        
        p <- p +
          # unadjusted (Wald) CI 
          ggplot2::geom_errorbar(
            orientation = "y",
            data = studyCIs,
            ggplot2::aes(y = y + y_offset, xmin = lower_unadj, xmax = upper_unadj),
            width     = diamond_height * 0.5,
            linewidth = 0.3,
            color     = "grey80",
            show.legend = FALSE
          ) +
          # raw reported estimate -> open circle 
          ggplot2::geom_point(
            data = studyCIs,
            ggplot2::aes(x = estimate_raw, y = y + y_offset),
            shape = 1L,
            color = "grey50"
          )
      }
      
      p <- p +
        ggplot2::geom_errorbar(
          orientation = "y",
          data = studyCIs,
          ggplot2::aes(
            y = y,
            xmin = lower,
            xmax = upper,
          ),
          width = diamond_height,
          show.legend = FALSE
        ) +
        # Adjusted (bias-corrected) point estimate -> filled circle
        ggplot2::geom_point(
          data = studyCIs,
          ggplot2::aes(x = estimate, y = y),
          shape = 16L
        )
    }     
      
    p <- p +
      ggplot2::geom_polygon(
        data = polygons[polygons$ci_exists == TRUE, ],
        ggplot2::aes(
          x = x,
          y = y,
          group = paste0(name, ".", id),
          fill = .data$color_hex
        ),
        show.legend = FALSE
      ) +
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
            labels = c(studyCIs$name, "", unique(polygons$name)),
            limits = c(min(spacing) - v_space/2, max(spacing) + v_space/2) # this fixes the problem where that the last polygon was always too close to the x axis! 
            
        ) +
        # ggplot2::labs(
        #     x = bquote(mu)
        # ) +
      ggplot2::scale_fill_identity()
    
    if (is.null(xlab)) {
        p <- p +
            ggplot2::labs(
                x = bquote(mu)
            )
    } else {
        p <- p +
            ggplot2::labs(
                x = xlab
            )
    }

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

    p
}


get_CI_new_methods <- function(
    cms,
    diamond_height,
    scale_diamonds
) {

    out <- lapply(seq_along(cms), function(r) {

        # Get  the cms
        cm <- cms[[r]]

        # Calculate polygons
        ci_exists <- if (all(is.na(cm$joint_cis))) FALSE else TRUE

        if (ci_exists) {
            # Calculate the polygons for the diamond
          
            # Basically gives back df 3 columns: x,y,id. If there is only one CI, then id == 1, then we have something like:
          #          row3(x3,+y3)
          #     row2 (x2,+y2)   row4(x4,+y4)
          # row1 (x1,0)                 row5(x5,0)
          #     row8(x2,-y2)   row6(x4,-y4)
          #          row7(x3,-y3)
          
          # since the polygon is bult from left to right following the row order
            polygons <- calculate_polygon(
                cm = cm,
                diamond_height = diamond_height,
                scale_diamonds = scale_diamonds
            )
            polygons$name <- cm$fun_name
            polygons$color <- r

        } else {
            polygons <- data.frame(
                x = NA_real_,
                y = NA_real_,
                id = NA_real_,
                name = cm$fun_name,
                color = r
            )
        }

        # Calculate the p-value at mu = 0
        p_0 <- data.frame(
            name = cm$fun_name,
            p_0 = cm$p_0[, 2L],
            stringsAsFactors = FALSE,
            row.names = NULL
        )

        p_max <- data.frame(
            name = cm$fun_name,
            mu = cm$p_max[, 1L],
            p_max = cm$p_max[, 2L],
            stringsAsFactors = FALSE,
            row.names = NULL
        )

        # Return
        list(
            CIs = polygons,
            p_0 = p_0,
            p_max = p_max
        )
    })
    
    
    # So out is a list containing each MA method used. e.g. out[[1]] = Edginton, out[[2]] = Fisher ... etc
    
    
    CIs <- do.call("rbind", lapply(out, "[[", i = 1L))
    p_0 <- do.call("rbind", lapply(out, "[[", i = 2L))
    p_max <- do.call("rbind", lapply(out, "[[", i = 3L))

    # Return a list of 3 objects 
    list(CIs = CIs, p_0 = p_0, p_max = p_max)
}


get_CI_old_methods <- function(
    cm,              
    diamond_height) 
  {


    # in previous version, there was a call to get_stats_others, and the CI comparison were recomputed. 
    #Now i just extract them 
  
   cis <- cm$comparison_cis
   p_0_vals <- cm$comparison_p_0
   
   
   #DIAMOND IS:
   # Left Corner: (Lower, 0)
   # Bottom Corner: (Est, -H/2) 
   # Right Corner: (Upper, 0)
   # Top Corner: (Est, +H/2)
   
  # get the estimates (center of CI)
   ests <- rowMeans(cis)
   
    # Create df
    t_ci <- t(cis) #row1 = contains all Lower bounds, row2 all lower bounds
    
    
    t_ci <- rbind(t_ci[1L, ], ests, t_ci[2L, ], ests) 
    
    # t_ci has:
    #Row 1: Lower bounds (Left corner x)
    #Row 2: Estimates (Bottom corner x)
    #Row 3: Upper bounds (Right corner x)
    #Row 4: Estimates (Top corner x)
    
    x <- c(t_ci) #flatten columnwise st we obtain:: Lower1, Est1, Upper1, Est1, Lower2, Est2, Upper2, Est2 ...
    
    #create y coordinates 
    y <- rep(
        c(0, -diamond_height / 2, 0, diamond_height / 2),
        times = length(ests)
    )
    
    method_names <- rownames(cis)

    CIs <- data.frame(
        x = x,
        y = y,
        id = 1L,
        name = rep(method_names, each = 4L),
        color = 0L,
        row.names = NULL,
        stringsAsFactors = FALSE
    )

    p_0 <- data.frame(
        name = method_names,
        p_0 = p_0_vals[, 2L],
        row.names = NULL,
        stringsAsFactors = FALSE
    )

    list(
        CIs = CIs,
        p_0 = p_0
    )
}

# Make a df that contains the data for the diamonds
calculate_polygon <- function(
    cm,
    diamond_height,
    scale_diamonds
) {

     if (all(is.na(cm$joint_cis))) return(NULL)
  
  
  
  no_gamma <- all(is.na(cm$gamma))
  
  # ------------------------------------------------------------------
  # Evaluate the p-value function at the key x-coordinates
  # (estimates, CI boundaries, maxima, minima) to get their heights.
  #
  # Point types:
  #   0 = lower CI boundary
  #   1 = local minimum (gamma)
  #   2 = estimate or local maximum
  #   3 = upper CI boundary
  # ------------------------------------------------------------------
  
  # p-value at each study estimate (evaluated at mu = estimate_i)
  f_estimates <- call_pfun(
    fun       = cm$p_fun,
    estimates = cm$estimates,
    SEs       = cm$SEs,
    mu        = cm$estimates,
    w         = cm$w,
    k         = cm$k_studies
  )
  
  nrep <- nrow(cm$joint_cis)
  
  pt_eval <- rbind(
    cbind(cm$p_max,                                              2L),  # maxima
    cbind(cm$estimates,     f_estimates,  rep(2L, length(cm$estimates))),  # estimates
    cbind(cm$joint_cis[, 1L], 0,          rep(0L, nrep)),                 # lower CI
    cbind(cm$joint_cis[, 2L], 0,          rep(3L, nrep))                  # upper CI
  )
  
  if (!no_gamma) {
    pt_eval <- rbind(pt_eval, cbind(cm$gamma, 1L))  # local minima
  }
  
  colnames(pt_eval) <- c("x", "y", "type")
  rownames(pt_eval) <- NULL
  
  # --- remove duplicates, keep only points inside the
  # CI, then sort by x so the polygon traces left to right.

  pt_eval <- pt_eval[!duplicated(pt_eval), ]
  
  inside_ci <- vapply(
    pt_eval[, "x"],
    function(xi) any(xi >= cm$joint_cis[, 1L] & xi <= cm$joint_cis[, 2L]),
    logical(1L)
  )
  pt_eval <- pt_eval[inside_ci, ]
  pt_eval <- pt_eval[order(pt_eval[, "x"]), ]
  
  
  # ---- scale y so the tallest point reaches diamond_height / 2.

  x    <- pt_eval[, "x"]
  y    <- pt_eval[, "y"]
  type <- pt_eval[, "type"]
  
  if (scale_diamonds) y <- y / max(y)
  y_scaled <- y * diamond_height / 2
  
  
  # ---  Assign a polygon id to each CI segment.---

  
  lower_idx <- which(type == 0L)
  upper_idx <- which(type == 3L)
  
  # in case we need to build more than one polygon (disconnected CI), assign to each one of those an id
  id <- integer(length(x))
  for (i in seq_along(lower_idx)) {
    id[lower_idx[i]:upper_idx[i]] <- i
  }
  
  # Each CI segment becomes a closed polygon: top edge (left to right)
  # then bottom edge (right to left, mirrored).

  do.call("rbind", lapply(unique(id), function(seg_id) {
    pts   <- id == seg_id
    x_seg <- x[pts]
    y_seg <- y_scaled[pts]
    n     <- length(x_seg)
    
    # Top edge: left to right at +y
    # Bottom edge: right to left at -y (excluding x_seg[-c(1L, n)] and the y since are the diamond edges and must be counted once)
    
    data.frame(
      x  = c(x_seg,               rev(x_seg[-c(1L, n)])),
      y  = c(y_seg,              -rev(y_seg[-c(1L, n)])),
      id = seg_id
    )
  }))
}





# ------------------------------------------------------------------------------
# Helper to call p_val function since we added weights and k
# ------------------------------------------------------------------------------

call_pfun <- function(fun, estimates, SEs, mu, w = NULL, k = NULL) {
  args <- names(formals(fun))
  
  # if 'w' and 'k' exist in the arguments and we have a valid w or k--> pass 
  has_w <- "w" %in% args && !is.null(w)
  has_k <- "k" %in% args && !is.null(k)
  
  if (has_w && has_k) {
    fun(estimates = estimates, SEs = SEs, mu = mu, w = w, k = k)
  } else if (has_w) {
    fun(estimates = estimates, SEs = SEs, mu = mu, w = w)
  } else if (has_k) {
    fun(estimates = estimates, SEs = SEs, mu = mu, k = k)
  } else {
    fun(estimates = estimates, SEs = SEs, mu = mu)
  }
}





#  ---------------------------------------------------------------------------- 
#   Add the bayesmeta object to the forest plot. 
# --------------------------------------------------------------------------------

add_bayes_forest <- function(p, bm, diamond_height, v_space, 
                             color = "#000000", label = "Bayesmeta",
                             mu_estimate = "median") {
  
  # Extract summary from the bayesmeta object
  s <- bm$summary
  mu_est <- s[mu_estimate, "mu"]
  lower  <- s["95% lower", "mu"]
  upper  <- s["95% upper", "mu"]
  
  # Get the current breaks and labels from the ggplot object's scale
  y_scale <- p$scales$get_scales("y")
  y_breaks <- y_scale$breaks
  y_labels <- y_scale$labels
  
  # Calculate new position (one step down using v_space)
  y_bayes <- min(y_breaks) - v_space
  height <- diamond_height / 2
  
  # Define polygon for the Bayesmeta diamond
  diamond <- data.frame(
    x = c(lower, mu_est, upper, mu_est),
    y = c(y_bayes, y_bayes + height, y_bayes, y_bayes - height)
  )
  
  # Add the diamond
  p <- p +
    ggplot2::geom_polygon(
      data = diamond,
      ggplot2::aes(x = x, y = y),
      fill = color,
      color = color
      )
  
  # Extend y-axis labels and breaks
  y_breaks <- c(y_breaks, y_bayes)
  y_labels <- c(y_labels, label)
  
  # Overwrite the scale
  suppressMessages({
    p <- p + ggplot2::scale_y_continuous(
      breaks = y_breaks,
      labels = y_labels,
      limits = c(y_bayes - v_space, max(y_breaks) + (v_space / 2))
      )
  })
  
  p
}
