library(patchwork)
library(devtools)
load_all()

n <- 20
mean <- 0
sd <- 1.5
shape <- 5
rate <- 5
heterogeneity <- c("none", "additive")
pValueFUN <- c("hMean", "Edgington", "Fisher")
distr <- "chisq"
level <- 0.95
diamond_height <- 0.5
v_space <- 1.5
studyNames <- NULL
xlim <- c(-5.4, 5)
scale_diamonds <- TRUE
rng_seed <- 42L
pValueFUN_args <- list(check_inputs = FALSE)

set.seed(rng_seed)
thetahat <- rnorm(n, mean = mean, sd = sd)
se <- rgamma(n, shape = shape, rate = rate)
pars <- list(
    thetahat = thetahat,
    se = se,
    het = heterogeneity,
    pValueFUN = pValueFUN,
    level = level,
    dist = distr,
    xlim = xlim,
    scale_diamonds = scale_diamonds
)

# Define a function that does the plots
make_plot <- function(pars) {
    # pval <- ggPvalueFunction(
    ggPvalueFunction(
        thetahat = pars$thetahat,
        se = pars$se,
        level = pars$level,
        heterogeneity = pars$het,
        distr = pars$distr,
        pValueFUN = pars$pValueFUN,
        pValueFUN_args = pValueFUN_args,
        xlim = pars$xlim
    )
    # forest <- ForestPlot(
    ForestPlot(
        thetahat = pars$thetahat,
        se = pars$se,
        level = pars$level,
        distr = pars$distr,
        pValueFUN = pars$pValueFUN,
        heterogeneity = pars$het,
        xlim = pars$xlim,
        scale_diamonds = pars$scale_diamonds
    )
    pval / forest
}

make_plot(pars)
