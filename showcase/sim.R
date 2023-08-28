library(patchwork)
library(devtools)
load_all()

n <- 9
mean <- 0
sd <- 2.5
shape <- 8
rate <- 5
heterogeneity <- c("none", "additive")
pValueFUN <- c("hMean", "Fisher")
distr <- "chisq"
level <- 0.95
diamond_height <- 0.5
v_space <- 1.5
studyNames <- NULL
xlim <- c(-6, 7)
scale_diamonds <- TRUE
v_space <- 1.5
diamond_height <- 0.5
show_studies <- TRUE
rng_seed <- 42L
pValueFUN_args <- list(check_inputs = FALSE)
drapery <- TRUE

set.seed(rng_seed)
thetahat <- rnorm(n, mean = mean, sd = sd)
se <- rgamma(n, shape = shape, rate = rate)

pars <- list(
    thetahat = thetahat,
    se = se,
    level = level,
    distr = distr,
    pValueFUN = pValueFUN,
    heterogeneity = heterogeneity,
    diamond_height = diamond_height,
    v_space = v_space,
    studyNames = studyNames,
    xlim = xlim,
    show_studies = show_studies,
    scale_diamonds = scale_diamonds,
    drapery = drapery
)

list2env(pars, globalenv())

# Define a function that does the plots
make_plot <- function(pars) {
    pval <- ggPvalueFunction(
    # ggPvalueFunction(
        thetahat = pars$thetahat,
        se = pars$se,
        level = pars$level,
        heterogeneity = pars$heterogeneity,
        distr = pars$distr,
        pValueFUN = pars$pValueFUN,
        pValueFUN_args = pValueFUN_args,
        xlim = pars$xlim
    )
    forest <- ForestPlot(
    # ForestPlot(
        thetahat = pars$thetahat,
        se = pars$se,
        level = pars$level,
        distr = pars$distr,
        pValueFUN = pars$pValueFUN,
        heterogeneity = pars$heterogeneity,
        diamond_height = pars$diamond_height,
        v_space = pars$v_space,
        studyNames = pars$studyNames,
        xlim = pars$xlim,
        show_studies = pars$show_studies,
        scale_diamonds = pars$scale_diamonds
    )
    pval$plot / forest$plot
}

make_plot(pars)
