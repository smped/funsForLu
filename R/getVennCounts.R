#' Prepare columns for plotting as a Venn Diagram
#'
#'@param data A data.frame or matrix
#'@param grpNames Optional. The names of the groups to be plotted
#'@param rmRexp Optional. A regular expression to be removed from the group names
#'
#'@return Returns a list with named objects ready to be plotted using \code{draw.quad.venn}, \code{draw.quintuple.venn} and associated functions
#'
#'@author Steve Pederson
#'
#'@import dplyr
#'@import VennDiagram
#'
#'@examples 
#'x <- data.frame(a = rbinom(20, 1, 0.5), b = rbinom(20, 1, 0.5))
#'getVennCounts(x)
#'
#'@rdname getVennCounts
#'@export
getVennCounts <- function(data, grpNames, rmRexp){
  
  stopifnot(any(is.matrix(data), is.data.frame(data))) # Must be 2 dimensional
  colClasses <- vapply(data, typeof, character(1))
  stopifnot(all(colClasses %in% c("integer", "logical"))) # No characters etc
  stopifnot(length(unique(colClasses)) == 1) # All the same data type
  if (colClasses[1] == "integer") stopifnot(all(max(data) == 1, min(data) == 0)) # No values that are non binary
  
  nCol <- ncol(data)
  nRow <- nrow(data)
  if(missing(grpNames)) {
    grpNames <- colnames(data)
    if(!missing(rmRexp)) grpNames <- gsub(rmRexp, "", grpNames)
  }
  else stopifnot(length(grpNames) == nCol)
  if (colClasses[1] != "logical") data <- apply(data, 2, as.logical)
  data <- as.data.frame(data)
  colnames(data) <- paste0("grp", 1:nCol)
  
  out <- lapply(data, sum)
  names(out) <- paste0("area", 1:nCol)
  out$n12 <- nrow(filter(data, grp1, grp2))
  if (nCol > 2){
    out$n13 <- nrow(filter(data, grp1, grp3))
    out$n23 <- nrow(filter(data, grp2, grp3))
    out$n123 <- nrow(filter(data, grp1, grp2, grp3))
  }
  if (nCol > 3){
    out$n14 <- nrow(filter(data, grp1, grp4))
    out$n24 <- nrow(filter(data, grp2, grp4))
    out$n34 <- nrow(filter(data, grp3, grp4))
    out$n124 <- nrow(filter(data, grp1, grp2, grp4))
    out$n134 <- nrow(filter(data, grp1, grp3, grp4))
    out$n234 <- nrow(filter(data, grp2, grp3, grp4))
    out$n1234 <- nrow(filter(data, grp1, grp2, grp3, grp4))
  }
  if (nCol > 4){
    out$n15 <- nrow(filter(data, grp1, grp5))
    out$n25 <- nrow(filter(data, grp2, grp5))
    out$n35 <- nrow(filter(data, grp3, grp5))
    out$n45 <- nrow(filter(data, grp4, grp5))
    out$n125 <- nrow(filter(data, grp1, grp2, grp5))
    out$n135 <- nrow(filter(data, grp1, grp3, grp5))
    out$n145 <- nrow(filter(data, grp1, grp4, grp5))
    out$n235 <- nrow(filter(data, grp2, grp3, grp5))
    out$n245 <- nrow(filter(data, grp2, grp4, grp5))
    out$n345 <- nrow(filter(data, grp3, grp4, grp5))
    out$n1235 <- nrow(filter(data, grp1, grp2, grp3, grp5))
    out$n1245 <- nrow(filter(data, grp1, grp2, grp4, grp5))
    out$n1345 <- nrow(filter(data, grp1, grp3, grp4, grp5))
    out$n2345 <- nrow(filter(data, grp2, grp3, grp4, grp5))
    out$n12345 <- nrow(filter(data, grp1, grp2, grp3, grp4, grp5))
  }
  
  out$category <- grpNames
  out
}