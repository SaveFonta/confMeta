#' @title Stouffer's method
#' @family p-value combination functions
#'
#' @description
#' Stouffer's method for combining \emph{p}-values across studies
#'
#' @inheritParams p_tippett
#' @param w Numeric vector of study weights. Defaults to \code{1/SEs} producing
#'     the same \emph{p}-value as Fixed Effects (Inverse Variance Weighted) meta-analysis
#' @param input_p Type of study-level \emph{p}-values used in the combination:
#'     \code{"greater"} (default) or \code{"less"}. Unlike the \emph{p}-value
#'     based methods, Stouffer's method does \strong{not} support
#'     \code{"two.sided"}, because mapping a \emph{z}-value to a \emph{p}-value
#'     and back requires a signed (directional) \emph{p}-value; a two-sided
#'     input discards the sign.  Supplying \code{input_p = "two.sided"} raises
#'     an error.
#'     
#'
#'
#' @details
#' Stouffer's \emph{z}-statistic for \eqn{n} independent studies is defined as
#' \deqn{z = \frac{\sum_{i=1}^n w_i z_i}{\sqrt{\sum_{i=1}^n w_i^2}},}
#' where \eqn{z_i} and \eqn{{w_i}} are individual study \emph{z}-values and
#' weights, respectively. Under the global null hypothesis, each \eqn{z_i} is
#' assumed to follow a standard normal distribution. Stouffer's \emph{z}-statistic
#' then also follows a standard normal distribution. The combined \emph{p}-value
#' is calculated as the probability of observing a value equal or greater than
#' \eqn{z} from this distribution: \deqn{p_S = \Pr(Z > z)}
#' 
#' With the default weights \eqn{w_i = 1/\sigma_i}, this is equivalent to a
#' fixed-effect inverse-variance meta-analysis.
#' 
#' @section Best-of-k adjustment:
#' Because Stouffer's method combines \emph{z}-values rather than
#' \emph{p}-values, the best-of-k adjustment is applied on the \emph{z}-scale.
#' For each study the reported \eqn{z_i} is mapped to its one-sided
#' \emph{p}-value, adjusted by \eqn{p \mapsto 1 - (1 - p)^{k_i}} (the
#' distribution function of the smallest of \eqn{k_i} independent uniforms),
#' and mapped back to an adjusted \emph{z}-value before pooling. The direction
#' of this mapping is determined by \code{input_p}. 
#'
#' Unlike the \emph{p}-value based methods, this requires a directional
#' \code{input_p} (\code{"greater"} or \code{"less"}): the sign of \eqn{z_i}
#' carries the direction of selection, which a two-sided \emph{p}-value would
#' discard.
#'
#' @inherit p_tippett return
#' 
#' 
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
#' p_stouffer(estimates, SEs, mu = 0, input_p = "greater", heterogeneity = "none")
#' 
#' Best-of-k adjustment: each study is the most significant of k = 3 experiments
#' p_stouffer(estimates, SEs, mu = 0, input_p = "greater", k = c(3, 3, 3))
p_stouffer <- function(
    estimates,
    SEs,
    mu = 0,
    heterogeneity = "none",
    phi = NULL,
    tau2 = NULL,
    check_inputs = TRUE,
    input_p = "greater",
    w = 1/SEs,
    output_p = "two.sided",
    k = rep(1, length(estimates))
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
      check_k_arg(k = k, estimates = estimates)
      check_w_arg(w = w, estimates = estimates)
    }

  
  # Stouffer cannot be two sided since the change from Z score to p value need to be one sided
  if (input_p == "two.sided") stop ("Stouffer supports only input_p 'greater' or 'less'")
  
  
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
    
  
    
    ## Apply best-of-k adjustment to Z-scores 
    # we need to map the z to p value based on input p, then apply adjustment, then map back 
    

    if (!all(k == 1L)) {
      z <- adjust_z_for_best_k(z = z, k = k, input_p = input_p)
    }
    

    ## compute weighted Stouffer's p-value
    zs <- colSums(w * z) / sqrt(sum(w^2))
    
    if (input_p == "less") {
      pstouffer <- stats::pnorm(zs, lower.tail = TRUE)
    } else {
      # Defaults to "greater" mapping
      pstouffer <- stats::pnorm(zs, lower.tail = FALSE)
    }
    
    
    if (output_p == "two.sided") {
      pstouffer <- 2 * pmin(pstouffer, 1 - pstouffer)
    }

    return(pstouffer)
}
