p_tippett <- function(
    thetahat,
    se,
    mu = 0,
    phi = NULL,
    tau2 = NULL,
    heterogeneity = "none",
    check_inputs = TRUE
) {

    # check inputs
    if (check_inputs) {
        check_inputs_p_value(
            thetahat = thetahat,
            se = se,
            mu = mu,
            heterogeneity = heterogeneity,
            phi = phi,
            tau2 = tau2
        )
    }

    # recycle `se` if needed
    if (length(se) == 1L) se <- rep(se, length(thetahat))

    # adjust se based on heterogeneity model
    se <- adjust_se(
      se = se,
      heterogeneity = heterogeneity,
      phi = phi,
      tau2 = tau2
    )

    # Get length
    n <- length(thetahat)

    # get the z-values
    z <- get_z(thetahat = thetahat, se = se, mu = mu)
    # convert them to p-values
    # p <- ReplicationSuccess::z2p(z, "two.sided")
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE) # faster than above
    # Calculate the Tippett statistic
    S_t <- apply(p, 2L, min)
    # Calculate the p-value using the beta distribution
    p <- stats::pbeta(q = S_t, shape1 = 1, shape2 = n, lower.tail = TRUE)
    # return
    return(p)
}
