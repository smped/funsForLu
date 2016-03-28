#' Get Confidence Intervals after the bootstrap procedure
#' 
#'@details Extract the confidence intervals from the provided samples across all tissues
#'
#'@param x The results from \code{\link{bootMeans}}
#'@param alpha Return (1 - alpha) confidence intervals
#'@param combined Combine all Tissues within the sample
#'
#'@return Returns a \code{tbl_df} with lower, median and upper quantiles for the interval. 
#'If \code{combined = TRUE} the \code{tbl_df} will only contain one row.
#'If \code{combined = FALSE} the \code{tbl_df} will contain one row per tissue.
#'
#'@author Steve Pederson
#'
#'@import dplyr
#'@import magrittr
#'
#'@rdname getCI
#'@export
getCI <- function(x, alpha, combined = FALSE){
  
  chkNames <- x %>% # Check that all tissues have the correct structure.
    vapply(function(x){
      "samples" %in% names(x)
    },
    logical(1))
  stopifnot(chkNames)
  
  if (!combined) {
    out <- x %>%
      lapply(extract2, "samples") %>%
      lapply(quantile, probs = c(alpha/2, 0.5, 1 - alpha/2)) %>%
      lapply(t) %>%
      lapply(as.data.frame) %>%
      bind_rows() %>%
      set_names(c("lower", "median", "upper")) %>%
      mutate(Tissue = names(x))
  } 
  else {
    out <- x %>%
      lapply(extract2, "samples") %>%
      unlist %>%
      quantile(probs = c(alpha/2, 0.5, 1 - alpha/2)) %>%
      as.list %>%
      as.data.frame %>%
      set_names(c("lower", "median", "upper"))
  }
  out
  
}