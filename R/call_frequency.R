#' Call Frequency
#'
#' \code{compute_cf} returns results from call frequency check
#'
#' @param x A matrix. Genotype names are stored in rows whereas marker names
#'  are stored in columns.
#' @param output Character vector with output options.
#' @param missing_value Character vector providing the encoding of missing elements.
#' @param call_threshold Numeric or complex vector specifying the call frequency
#'  threshold.
#' @return If \code{output} is "marker_names" a character vector with marker
#'  names that have passed the quality check will be returned. If \code{output}
#'  is "marker_callfreq" a numeric vector with the call frequency at each
#'  marker locus will be returned.
#' @examples
#'  # Load a matrix with SNP genotypes encoded as numeric values
#'  data(marker_numeric)
#'
#'  # Compute the call frequency
#'  compute_cf(marker_numeric, output = "marker_callfreq",
#'             missing_value = NA_real_)
#'
#'  # Retrieve the marker names of all loci with a call frequency >= 0.9.
#'  compute_cf(marker_numeric, output = "marker_names",
#'             missing_value = NA_real_, call_threshold = 0.9)
#' @export
compute_cf <- function(x, output = c("marker_names", "marker_callfreq"),
                       missing_value = "??", call_threshold = NULL) {

  if (storage.mode(x) != typeof(missing_value)) {
    stop("Missing values do not have the same class as genotypes")
  }
  if (output == "marker_names" && is.null(call_threshold)) {
    stop("'call_threshold' is undefined")
  }
  # Compute the call frequency for each marker.
  if (storage.mode(x) == "double") {
    callfreq <- 1 - colMeans(is.na(x))
  } else if (storage.mode(x) == "character") {
    callfreq <- 1 - colMeans(x == missing_value)
  }
  output <- match.arg(output)
  switch(EXPR = output,
         marker_names = colnames(x)[callfreq >= call_threshold],
         marker_callfreq = callfreq)
}
