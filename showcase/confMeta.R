library(devtools)
load_all()

n <- 5
mean <- 0
sd <- 2.5
shape <- 8
rate <- 5
# heterogeneity <- c("none", "additive")
# pValueFUN <- c("hMean", "Fisher")
# distr <- "chisq"
level <- 0.95
study_names <- NULL
# xlim <- c(-6, 7)
# scale_diamonds <- TRUE
# v_space <- 1.5
# diamond_height <- 0.5
# show_studies <- TRUE
rng_seed <- 42L
# pValueFUN_args <- list(check_inputs = FALSE)
# drapery <- TRUE

set.seed(rng_seed)
estimates <- rnorm(n, mean = mean, sd = sd)
SEs <- rgamma(n, shape = shape, rate = rate)
