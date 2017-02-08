#' Heterozygosity
#'
#' \code{compute_het} returns results from heterozygosity check
#'
#' @param x A matrix. Genotype names are stored in rows whereas marker names
#'  are stored in columns. No missing values are allowed.
#' @param output Character vector with output options.
#' @param hetThresh Numeric or complex vector specifying the heterozygosity
#'  threshold.
#' @param na_coding Scalar providing the encoding of missing elements.
#' @return If \code{output} is "marker_names" a character vector with marker
#'  names that have passed the quality check will be returned. If \code{output}
#'  is "marker_heterozygosity" a numeric vector with the heterozygosity at each
#'  marker locus will be returned.
#' @examples
#'  # Load a matrix with SNP genotypes encoded as numeric values
#'  data(imp_snps)
#'
#'  # Return the names of all markers with a heterozygosity of less than 10%.
#'  compute_het(imp_snps, output = "marker_names", het_threshold = 0.1)
#'
#'  # Return the heterozygosity at each locus.
#'  compute_het(imp_snps, output = "marker_heterozygosity", het_threshold = 0)
#' @export
### Remove markers with high heterozygosity levels.
compute_het <- function(x, output = c("marker_names", "marker_heterozygosity"),
                        het_threshold = 0, na_coding = "??") {

  # Recode missing values that are actually encoded as 'NA' to avoid conflicts.
  if (anyNA(x)) {
    if (storage.mode(x) == "double") {
      na_coding <- 99
    } else if (storage.mode(x) == "character") {
      na_coding <- "??"
    }
    x[is.na(x)] <- na_coding
  }

  if (is.null(colnames(x))) stop("Assign colnames to x")
  output <- match.arg(output)
  allgeno <- unique(c(x[x != na_coding]))

  if (storage.mode(x) == "double") {
    put <- sort(allgeno)
    if (mean(put[c(1, 3)]) != put[2]) {
        stop("Value for heterozygotes must be equal to the mean of homozygotes")
      } else {
        storage.mode(x) <- "character"
        storage.mode(put) <- "character"
        x[x == put[1]] <- "AA"
        x[x == put[2]] <- "AB"
        x[x == put[3]] <- "BB"
      }
  }
  allgeno <- unique(c(x))

  if (!all(nchar(allgeno) == 2)) {
    stop(cat("Alleles must be encoded by two letters."))
  }
  # Split genotypes into their alleles...
  namesplit <- strsplit(allgeno, split = "")
  # ... and extract heterozygous genotypes.
  het <- allgeno[unlist(lapply(namesplit, function(x) x[1] != x[2]))]
  # Set all missing values to 'NA_character_'.
  x[x %in% na_coding] <- NA_character_
  # Code heterozygous genotypes as '1' and all others as '0'.
  x[x %in% het] <- "1"
  x[x != "1"] <- "0"
  snp_het <- colSums(x == "1", na.rm = TRUE) / nrow(x)
  switch(EXPR = output,
         marker_heterozygosity = snp_het,
         marker_names = colnames(x[, snp_het <= het_threshold]))
}
