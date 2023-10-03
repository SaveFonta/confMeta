library(patchwork)
library(devtools)
load_all()

n <- 9
mean <- 0
sd <- 1.5
shape <- 8
rate <- 5
heterogeneity <- c("none", "additive")
p_fun <- p_edgington
# pValueFUN <- c("hMean", "Fisher")
distr <- "chisq"
conf_level <- 0.95
diamond_height <- 0.5
v_space <- 1.5
study_names <- NULL
xlim <- c(-6, 7)
scale_diamonds <- TRUE
show_studies <- TRUE
rng_seed <- 42L
pValueFUN_args <- list(check_inputs = FALSE)
drapery <- TRUE

set.seed(rng_seed)
estimates <- rnorm(n, mean = mean, sd = sd)
estimates1 <- rnorm(n, mean = mean, sd = sd)
estimates2 <- rnorm(n, mean = mean, sd = sd)
SEs <- rgamma(n, shape = shape, rate = rate)
SEs1 <- rgamma(n, shape = shape, rate = rate)
SEs2 <- rgamma(n, shape = shape, rate = rate)
fun <- p_hmean
fun1 <- p_edgington
fun2 <- p_fisher
ell <- list(approx = FALSE, rando = "hi")
study_names <- NULL

type <- c("p", "forest")

cm <- confMeta(
    estimates = estimates,
    SEs = SEs,
    conf_level = conf_level,
    fun = p_hmean
)

cm1 <- confMeta(
    estimates = estimates,
    SEs = SEs,
    conf_level = conf_level,
    fun = p_edgington
)

cm2 <- confMeta(
    estimates = estimates,
    SEs = SEs,
    conf_level = conf_level,
    fun = p_fisher
)

cms <- list(cm, cm1, cm2)

ggplot2::autoplot(cm)

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
