#' @title Concatenate the column names from a logical data.frame
#' 
#' @description This checks that a data.frame is all logical values, 
#' then concatenates the column names using the given separator.
#' 
#' @param df The \code{data.frame}
#' @param collapse The separator to use between the character variables in the column names
#' 
#' @details This is useful for turning sets of logical vectors into concatentaed strings 
#' able to be used for assigning values to bins.
#' 
#' @return A character vector
#' 
#' @examples 
#' t <- c(TRUE, FALSE)
#' x <- data.frame(Type1 = sample(t, 10, TRUE), Type2 = sample(t, 10, TRUE))
#' binsFromCols(x)
#' 
#' @author Steve Pederson
#' 
#' @export
binsFromCols <- function(df, collapse="_"){
  
  # Ensure the input is a df with only logical vectors
  stopifnot(vapply(df, is.logical, logical(1)))
  
  grps <- names(df)
  out <- apply(df, 1, function(x){paste0(grps[x], collapse = collapse )})
  out[out == ""] <-"none"
  out
  
}