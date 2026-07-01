#' @param k A numeric vector specifying the number of candidate studies 
#'     from which each reported estimate was selected as the "best" (most significant). 
#'     Defaults to `1` for all studies, which implies no selection bias (standard 
#'     unadjusted meta-analysis). Must be of the same length as `estimates` and 
#'     `SEs` (or length 1, which will be recycled).
