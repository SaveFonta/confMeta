#' @title Stouffer's method
#' @family p-value combination functions
#'
#' @description
#' Stouffer's method for combining \emph{p}-values across studies
#'
#' @inheritParams p_tippett
#' @param w Numeric vector of study weights. Defaults to \code{1/SEs} producing
#'     the same *p*-value as meta-analsis
#'
#'
#' @details
#' Stouffer's *z*-statistic for \eqn{k} studies is defined as
#' \deqn{z = \frac{\sum_{i=1}^k w_i z_i}{\sqrt{\sum_{i=1}^k w_i^2}},}
#' where \eqn{i_i} and \eqn{{w_i}} are individual study \emph{z}-values and
#' weights, respectively. Under the global null hypothesis, each \eqn{z_i} is
#' assumed to follow a standard normal distribution. Stouffer's *z*-statistic
#' then also follows a standard normal distribution. The combined \emph{p}-value
#' is calculated as the probability of observing a value equal or greater than
#' \eqn{z} from this distribution: \deqn{p_S = \Pr(Z > z)}
#'
#' @inherit p_tippett return
#'
#' @importFrom stats pnorm
#' @export
#'
#'
#'
#' @references
#' Stouffer SA, et al. *The American Soldier*. Princeton University Press, 1949.
#'  \doi{10.2307/2572105} 
#'   
#' Held, L, Hofmann, F, Pawel, S. (2025). A comparison of combined *p*-value
#' functions for meta-analysis. *Research Synthesis Methods*, 16:758-785.
#' \doi{10.1017/rsm.2025.26}
#'  
#' @examples
#' estimates <- c(0.1, 0.2, 0.3)
#' SEs <- c(0.05, 0.05, 0.1)
#' p_stouffer(estimates, SEs, mu = 0, heterogeneity = "none")
p_stouffer <- function(
    estimates,
    SEs,
    mu = 0,
    heterogeneity = "none",
    phi = NULL,
    tau2 = NULL,
    check_inputs = TRUE,
    w = 1/SEs,
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

    ## get the z-values
    z <- get_z(estimates = estimates, SEs = SEs, mu = mu)

    ## compute weighted Stouffer's p-value
    zs <- colSums(w * z) / sqrt(sum(w^2))
    pstouffer <- stats::pnorm(zs, lower.tail = FALSE)

    if (output_p == "two.sided") {
      pstouffer <- 2 * pmin(pstouffer, 1 - pstouffer)
    }

}
