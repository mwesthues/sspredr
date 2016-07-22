#' Call Frequency
#'
#' \code{compute_cf} returns results from call frequency check
#'
#' @param x A matrix. Genotype names are stored in rows whereas marker names
#'  are stored in columns.
#' @param output Character vector with output options.
#' @param missing Character vector providing the encoding of missing elements.
#' @param callThresh Numeric or complex vector specifying the call frequency
#'  threshold.
#' @return If \code{output} is "markerNames" a character vector with marker
#'  names that have passed the quality check will be returned. If \code{output}
#'  is "markerCallfreq" a numeric vector with the call frequency at each
#'  marker locus will be returned.
#' @examples
#'  # Load a matrix with SNP genotypes encoded as numeric values
#'  data(marker_numeric)
#'
#'  # Compute the call frequency
#'  compute_cf(marker_numeric, output = "markerCallfreq", missing = NA_real_)
#'
#'  # Retrieve the marker names of all loci with a call frequency >= 0.9.
#'  compute_cf(marker_numeric, output = "markerNames", missing = NA_real_,
#'             callThresh = 0.9)
#' @export
compute_cf <- function(x, output = c("markerNames", "markerCallfreq"),
                       missing = "??", callThresh = NULL) {
  # class(x) = matrix, rows = genotypes, columns = marker loci
  # Names of all genotypes present in the marker data matrix.
  # output:   a) markerName = vector with markers who passed the criterion
  #           b) markerCallfreq = named vector with call frequencies for
  #               each vector
  # callThresh: call frequency threshold
  output <- match.arg(output)

  if (storage.mode(x) != typeof(missing)) {
    stop("Missing values do not have the same class as genotypes")
  }
  if (output == "markerNames" && is.null(callThresh)) {
    stop("'callThresh' is undefined")
  }
  # Compute the call frequency for each marker.
  if (storage.mode(x) == "double") {
    callfreq <- 1 - colMeans(is.na(x))
  } else if (storage.mode(x) == "character") {
    callfreq <- 1 - colMeans(x == missing)
  }
  switch(EXPR = output,
         markerNames = colnames(x)[callfreq >= callThresh],
         markerCallfreq = callfreq)
}
