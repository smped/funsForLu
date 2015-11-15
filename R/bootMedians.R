#' @title  Bootstrap two gene sets and compare their values
#'
#' @description A function for comparing two sets of genes without relying on any distributional assumptions.
#'
#' @param data A data frame containing all the data required
#' @param testIds A \code{character vector} with test set of Ids
#' @param refIds A \code{character vector} with the reference setof Ids
#' @param idCol The column in \code{data} containing the Ids in the vectors \code{testIds} and \code{refIds}. Can be specified as an integer position or as a character (regular expression).
#' @param binCol The column in \code{data} containing the bin allocations for each gene. Can also be specified as an integer or by name.
#' @param valCols Regular expression or integers defining the columns in \code{data} containing the values of interest.
#' @param nGenes \code{integer}. The number of genes to sample at each iteration. Values greater than the number of testIds will automatically be capped at the number of testIds
#' @param nBoot \code{integer}. The number of bootstrap iterations to be performed
#' @param minGenes \code{integer}. The minimum number of IDs required to conduct a bootstrap procedure with any meaning.
#' @param ... Passed to the function \code{mean} internally
#' @param na.rm \code{logical}. Also passed internally to the function \code{mean}
#' @param replace \code{logical}. Should the bootstrap use sampling with replacement (\code{replace = TRUE}) or without
#' @param maxP The maximum probability (weight) allowed for an individual gene in the reference set. 
#' Defaults to \code{1/nGenes}
#'
#' @details
#'This function breaks the supplied \code{data.frame} into two sets of test IDs & reference IDs.
#'The \code{data.frame} must contain a column (\code{binCol}) which classifies each ID into a bin.
#'The probabilities of bin membership in the test IDs are then used for sampling during the bootstrap procedure.
#'
#'The values to be bootstrapped must be specified in the argument \code{valCols}, 
#'and this can be a regular expression used to extract a set of columns, 
#'or integers specifying exact columns in the supplied \code{data.frame}.
#'
#'The function itself will sample the same number of IDs (\code{nGenes}) from each dataset, 
#'based on the probabilities of bin membership in the test dataset.
#'At each bootstrap iteration, the median values for each column specified will be returned from both datasets,
#'with the reference values then subtracted from the tested values.
#'This allows direct comparison of these values as they will be drawn from similar
#'distributions based on the binning variable used.
#'
#'If any genes have a probability of being resampled > \code{maxP} they may exert undue influence on the results.
#'If any are found the process will stop to allow removal of this grouping.
#'Alternatively, the value for \code{maxP} can be reset up to a maximum of 1, 
#'which would represent maximum permissability.
#'
#' @author Steve Pederson
#'
#' @return 
#' A \code{list} with components:
#' \itemize{
#'  \item \code{$samples} The sampled differences in the median values
#'  \item \code{$p} The proportion of sampled differences which are > 0
#'  \item \code{$nGenes} The number of genes sampled at each bootstrap iteration
#'  \item \code{$nBoot} The number of bootstrap iterations
#'  \item \code{sampleSizes} A named numeric vector with the sample sizes for each dataset
#'  \item \code{$testBins} The distributions of genes amongst the binning variable in the set of test IDs
#'  \item \code{$refBins} The distributions of genes amongst the binning variable in the set of reference IDs.
#'  The final column represents the sampling probability for each individual gene in the corresponding bin
#'}
#' @import dplyr
#'
#' @export
bootMedians <- function(data, testIds, refIds, 
                        idCol = 1L, binCol = "lengthBin", valCols = "TPM",
                        nGenes = 1000L, nBoot = 100L, 
                        minGenes = 200L, ..., na.rm=TRUE, replace = TRUE, maxP){
  
  stopifnot(is.data.frame(data))
  stopifnot(is.character(testIds))
  stopifnot(is.character(refIds))
  stopifnot(is.numeric(nGenes))  
  stopifnot(is.numeric(nBoot))
  stopifnot(is.logical(na.rm), is.logical(replace))
  
  nGenes <- as.integer(nGenes)
  nBoot <-as.integer(nBoot)
  if (length(nBoot) > 1){
    nBoot <- nBoot[1]
    message("length(nBoot) > 1. Only the first value will be used.")
  }
  if (nBoot > 1e6) stop("You cannot bootstrap more than a million times. Sorry...")
  
  stopifnot(any(is.numeric(idCol), is.character(idCol)))
  if (length(idCol) > 1) {
    message("length(idCol) > 1. Only the first value will be used.")
    idCol <- idCol[1]
  }
  if (is.numeric(idCol))  { 
    stopifnot(idCol <= ncol(data))
    idCol <- as.integer(idCol)
  }
  if (is.character(idCol)){
    idCol <- grep(idCol, colnames(data), value = TRUE)
    if (length(idCol) != 1) stop("Unable to correctly identify idCol.\n Does it match multiple or no columns?")
  }
  data <- mutate(data, ID = data[[idCol]])   
  
  stopifnot(any(is.numeric(binCol), is.character(binCol)))
  if (length(binCol) > 1) {
    message("length(binCol) > 1. Only the first arguments will be used")
    binCol <- binCol[1]
  }
  if (is.numeric(binCol))  { 
    stopifnot(binCol <= ncol(data))
    binCol <- as.integer(binCol)
  }
  if (is.character(binCol)){
    binCol <- grep(binCol, colnames(data), value = TRUE)
    if (length(binCol) != 1) stop("Unable to correctly identify binCol.\n Does it match multiple or no columns?")
  }
  data <- mutate(data, bins = as.factor(data[[binCol]]))
  
  stopifnot(any(is.character(valCols), is.numeric(valCols)))
  if (is.numeric(valCols)){
    valCols <- abs(as.integer(valCols))
    stopifnot(valCols < ncol(data))
    valCols <- colnames(data)[valCols]
  }
  if (is.character(valCols)){
    valCols <- grep(valCols, colnames(data), value = TRUE)
  }
  
  idMatch <- c(test = sum(data$ID %in% testIds),
               ref  = sum(data$ID%in% refIds))
  if (nGenes > min(idMatch)) {
    # If the number of genes to be sampled from (n) is less than that requested (m), 
    # reset the value to m = n - 3
    nGenes <- min(idMatch) - 3
    message("NB: The maximum number of genes able to be bootstrapped is ", nGenes, ".\n")
    if (nGenes < minGenes) stop("Your test dataset is < ", minGenes, " IDs.\n",
                                "If this is expected behaviour, ",
                                "rerun the function changing the parameter minGenes to be below this value.\n",
                                "The number of matching IDs in your test & reference sets are ", paste(idMatch, collapse = " & "), 
                                " respectively")
  }
  
  if(missing(maxP)) {
    maxP <- 1 / nGenes
  }
  else {
    stopifnot(is.numeric(maxP), maxP <= 1) # Ensure a numeric value <= 1
  }
  
  # Now just get the key information from the supplied data
  data <- select(data, ID, bins, one_of(valCols))
  
  # Set up the test data frame
  testData <- filter(data, ID %in% testIds)
  testBins <- summarise(group_by(testData, bins), count = n())
  testBins <- mutate(ungroup(testBins), p = count / sum(count))
  
  # Set up the reference data frame
  refData <- filter(data, ID %in% refIds)
  refBins <- summarise(group_by(refData, bins), count = n())
  refBins <- mutate(ungroup(refBins), 
                    p = count / sum(count), 
                    effP = testBins$p / count) # This is the probability/gene
  
  # Assign bin weights for the reference dataset to ensure the matching sampling distributions
  wRef <- refBins$effP[match(refData$bins, refBins$bin)]
  refBins <- mutate(refBins, effP = effP / sum(wRef)) # Rescale effP to ensure a sum to 1
  
  # Check for any points with potentially excessive influence.
  # This defaults to checking for points with a probability of inclusion > 1 for each iteration
  if (max(refBins$effP) > maxP) stop("Genes from one or more groups have probabilities of being resampled >",
                                     maxP, "\n", refBins)
  
  # The sampling function
  sampleBoth <- function(testData, refData, nGenes, testBins, replace, wRef){
    
    # Take the samples
    testSamp <- sample_n(testData, nGenes, replace)
    refSamp <- sample_n(refData, nGenes, replace, weight = wRef)
    # Find the median expression value & find the difference
    testMeds <- summarise_each(testSamp, funs(median), one_of(valCols))
    refMeds <- summarise_each(refSamp, funs(median), one_of(valCols))
    meds <- testMeds - refMeds
    meds
    
  }
  
  # Now run the tilted bootstrap
  out <- replicate(nBoot, 
                   sampleBoth(testData, refData, nGenes, testBins, replace, wRef), 
                   simplify = FALSE)
  out <- bind_rows(out)
  p <- vapply(out,function(x){mean(x > 0, na.rm =na.rm, ...)}, numeric(1))
  
  # And the results
  return(list(samples = out,
              p = p,
              nGenes = nGenes,
              nBoot = nBoot,
              sampleSizes = c(test = length(testIds), 
                              ref = length(refIds)),
              testBins = testBins,
              refBins = refBins))
  
}

