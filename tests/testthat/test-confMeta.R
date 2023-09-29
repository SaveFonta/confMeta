test_that("argument checks work properly.", {

    n <- 5
    e <- estimates <- rnorm(n)
    s <- SEs <- rgamma(n, 5, 5)
    c <- conf_level <- 0.95
    fun <- p_edgington
    f <- fo <- formals(fun)

    # should return an object
    cl1 <- quote({
        confMeta(
            estimates = estimates,
            SEs = SEs,
            conf_level = conf_level,
            fun = fun
        )
    })
    expect_true(is.object(eval(cl1)))

    # =========================================================================
    # Try wrong inputs when creating a confMeta object
    # =========================================================================

    # estimates is character
    estimates <- letters[1:n]
    expect_error(eval(cl1))
    estimates <- e

    # SEs is character
    SEs <- letters[1:5]
    expect_error(eval(cl1))
    SEs <- s

    # NAs in estimates/SEs
    estimates[1] <- NA
    expect_error(eval(cl1))
    estimates[1] <- NaN
    expect_error(eval(cl1))
    estimates <- e
    SEs[1] <- NA
    expect_error(eval(cl1))
    SEs <- s

    # confidence level has wrong type
    conf_level <- "a"
    expect_error(eval(cl1))
    conf_level <- c

    # confidence level has wrong length
    conf_level <- rep(conf_level, 2)
    expect_error(eval(cl1))
    conf_level <- c

    # confidence level has invalid value
    conf_level <- 1.5
    expect_error(eval(cl1))
    conf_level <- c

    # confidence level is not finite
    conf_level <- NA_real_
    expect_error(eval(cl1))
    conf_level <- c

    # length mismatch
    estimates <- estimates[-1]
    expect_error(eval(cl1))
    estimates <- e

    # function with one of the required arguments missing
    formals(fun) <- formals(fun)[-1]
    expect_error(eval(cl1))
    formals(fun) <- f

    # Remove the default value for one of the non-required arguments
    fo["approx"] <- NULL
    fo <- append(fo, alist(approx = ))
    formals(fun) <- fo
    expect_error(eval(cl1))
    formals(fun) <- f

    # =========================================================================
    # Try the validate_confMeta function on some wrong fields
    # =========================================================================

    c <- cm <- eval(cl1)
    cl2 <- quote(validate_confMeta(cm))

    # fields missing
    cm[c("estimates", "SEs")] <- NULL
    expect_error(eval(cl2))
    cm <- c

    # necessary function arguments missing
    formals(cm$p_fun)[c("estimates", "mu")] <- NULL
    expect_error(eval(cl2))
    cm <- c

    # additional function arguments have no default
    fo["approx"] <- NULL
    formals(cm$p_fun) <- append(fo, alist(approx = ))
    expect_error(eval(cl2))
    fo <- f
    cm <- c

    # NAs/NaNs in estimates / SEs / conf_level
    cm$conf_level <- NA_real_
    expect_error(eval(cl2))
    cm <- c
    cm$estimates[1] <- NA_real_
    expect_error(eval(cl2))
    cm <- c
    cm$SEs[1] <- NA_real_
    expect_error(eval(cl2))
    cm <- c

})
