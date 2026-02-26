#' @title Pearson's method
#' @family p-value combination functions
#'
#' @description
#' Pearson's method for combining \emph{p}-values across studies. 
#' The method transforms the complement of the individual study \emph{p}-values 
#' using the natural logarithm and evaluates the sum against a chi-squared 
#' distribution.
#'
#' @inheritParams p_tippett
#'
#' @details
#' The Pearson test statistic for \eqn{k} studies is defined as:
#' \deqn{g = -2 \sum_{i=1}^k \log(1 - p_i)}
#' 
#' Under the global null hypothesis, each \eqn{p_i} is assumed to follow a 
#' \eqn{Unif(0, 1)} distribution. The test statistic \eqn{g} therefore follows a 
#' chi-squared distribution with \eqn{2k} degrees of freedom: \eqn{\chi^2_{2k}}. 
#' The combined \emph{p}-value, \eqn{p_P}, is calculated as the probability of observing a 
#' value less than or equal to \eqn{g} from this distribution:
#' \deqn{p_P = P(\chi^2_{2k} \leq g)}
#'
#' **Important note on orientation:** Unlike Edgington's method, Pearson's method 
#' is *not* orientation-invariant. The combined \emph{p}-value depends on 
#' the direction of the one-sided \emph{p}-values (controlled by the 
#' \code{input_p} argument). 
#' 
#' Specifically, Pearson's and Fisher's methods are mirrored. Computing 
#' the Pearson combined \emph{p}-value for the "greater" alternative is equal to 
#' 1 minus the Fisher combined \emph{p}-value for the "less" alternative.
#'
#' @inheritSection p_tippett Output p-value
#'
#' @inherit p_tippett return
#'
#' @importFrom stats pnorm pchisq
#'
#' @export
#'
#' @references
#' Pearson K. On a method of determining whether a sample of size n supposed to 
#' have been drawn from a parent population having a known probability integral 
#' has probably been drawn at random. *Biometrika*, 25:379-410, 1933.
#' 
#' Held L, Hofmann F, Pawel S. A comparison of combined p-value functions for 
#' meta-analysis. *Research Synthesis Methods*, 16:758-785, 2025.
#'
#' @examples
#' # Simulating estimates and standard errors
#' n <- 15
#' estimates <- rnorm(n)
#' SEs <- rgamma(n, 5, 5)
#'
#' # Set up a vector of means under the null hypothesis
#' mu <- seq(
#'    min(estimates) - 0.5 * max(SEs),
#'    max(estimates) + 0.5 * max(SEs),
#'    length.out = 100
#' )
#'
#' # Using Pearson's method to calculate the combined p-value
#' p_pearson(
#'      estimates = estimates,
#'      SEs = SEs,
#'      mu = mu,
#'      heterogeneity = "none",
#'      output_p = "two.sided",
#'      input_p = "greater"
#' )
#' 

p_pearson <- function(
    estimates,
    SEs,
    mu = 0,
    heterogeneity = "none",
    phi = NULL,
    tau2 = NULL,
    check_inputs = TRUE,
    input_p = "greater",
    output_p = "two.sided"
    ) {
  
    # check inputs
    if (check_inputs) {
        check_inputs_p_value(
            estimates = estimates,
            SEs = SEs,
            mu = mu,
            heterogeneity = heterogeneity,
            phi = phi,
            tau2 = tau2
        )
      check_output_p_arg(output_p = output_p)
    }

    # recycle `se` if needed
    if (length(SEs) == 1L) SEs <- rep(SEs, length(estimates))

    # adjust se based on heterogeneity model
    SEs <- adjust_se(
        SEs = SEs,
        heterogeneity = heterogeneity,
        phi = phi,
        tau2 = tau2
    )

    # Get lengths
    n <- length(estimates)

    # implement alternatives
    z <- get_z(estimates = estimates, SEs = SEs, mu = mu)
    
    p <- switch(input_p,
                "two.sided" = 2 * stats::pnorm(abs(z), lower.tail = FALSE),
                "greater"   = stats::pnorm(z, lower.tail = FALSE),
                "less"      = stats::pnorm(z, lower.tail = TRUE),
                stop("input_p must be 'greater','less','two.sided'")
    )
    
    p <- as.matrix(p)
    
    # Pearson exact calculation: P(ChiSq_2k <= -2 * sum(log(1-p)))
    ppearson <- stats::pchisq(
        q = -2 * colSums(log(1 - p)),
        df = 2 * n,
        lower.tail = TRUE
    )
    if (output_p == "two.sided" && input_p != "two.sided") {
      ppearson <- 2 * pmin(ppearson, 1 - ppearson)
    }
    
    return(ppearson)
}
