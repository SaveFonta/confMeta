library(patchwork)

n <- 15
mean <- 0
sd <- 1.3
shape <- 5
rate <- 5
heterogeneity <- c("none", "additive")
pValueFUN <- c("hMean", "Pearson", "Edgington", "Fisher")
distr <- "chisq"
level <- 0.95
diamond_height <- 0.5
v_space <- 1.5
studyNames <- NULL
xlim <- c(-5, 5)
rng_seed <- 43L
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
    xlim = xlim
)

# Define a function that does the plots
make_plot <- function(pars) {
    pval <- ggPvalueFunction(
        thetahat = pars$thetahat,
        se = pars$se,
        level = pars$level,
        heterogeneity = pars$het,
        distr = pars$distr,
        pValueFUN = pars$pValueFUN,
        pValueFUN_args = pValueFUN_args,
        xlim = pars$xlim

    )
    forest <- ForestPlot(
        thetahat = pars$thetahat,
        se = pars$se,
        level = pars$level,
        distr = pars$distr,
        pValueFUN = pars$pValueFUN,
        heterogeneity = pars$het,
        xlim = pars$xlim,
    )
    pval / forest
}

make_plot(pars)
