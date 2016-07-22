#' Minor Allele Frequency
#'
#' \code{compute_maf} returns results from minor allele frequency check
#'
#' @param x A matrix. Genotype names are stored in rows whereas marker names
#'  are stored in columns.
#' @param output Character vector with output options.
#' @param missing Character vector providing the encoding of missing elements.
#' @param mafThresh Numeric or complex vector specifying the minor allele
#'  frequency threshold.
#' @return If \code{output} is "markerNames" a character vector with marker
#'  names that have passed the quality check will be returned. If \code{output}
#'  is "markerMAF" a numeric vector with the minor allele frequency at each
#'  marker locus will be returned. If \code{output} is "genoList" a list of two
#'  elements will be returned, which provides the minor allele at each locus as
#'  well as the major allele. The speed of computations is considerably higher
#'  if the elements of \code{x} are supplied as character values since an
#'  internal class conversion will be applied otherwise.
#' @examples
#'  # Load a matrix with SNP genotypes encoded as numeric values
#'  data(marker_numeric)
#'
#'  # Return the names of all polymorphic markers.
#'  compute_maf(marker_numeric, output = "markerNames", missing = NA_real_,
#'              mafThresh = 0)
#'
#'  # Set one locus to monomorphic
#'  marker_mono <- marker_numeric
#'  marker_mono[, 1] <- 1
#'  compute_maf(marker_numeric, output = "markerNames", missing = NA_real_,
#'              mafThresh = 0)
#'
#'  # Load a matrix with SNP genotypes encoded as character values
#'  data(marker_character)
#'
#'  # Return the minor allele frequency at each locus with MAF >= 0.05.
#'  compute_maf(marker_character, output = "markerMAF", missing = "??",
#'              mafThresh = 0.05)
#'
#'  # Return the genotype of the major as well as of the minor allele at each
#'  # locus.
#'  compute_maf(marker_character, output = "genoList", missing = "??",
#'              mafThresh = 0)
#' @export
compute_maf <- function(x, output = c("markerNames", "markerMAF", "genoList"),
                        missing = "??", mafThresh = 0) {
  if (class(x) != "matrix") {
    stop("Input has to be of class 'matrix'")
  }
  if (class(missing) != class(c(x))) {
    stop("Classes of 'missing' and 'x' must be the same.")
  }
  # Recode SNP markers to class 'character' if they are not yet incoded as such.
  if (class(c(x)) != "character") {
    put <- sort(na.exclude(unique(c(x)[!c(x) %in% missing])))
      storage.mode(x) <- "character"
      x[x %in% missing] <- "??"
      missing <- "??"
    if (length(put) == 2) {
      x[x == put[1]] <- "AA"
      x[x == put[2]] <- "BB"
    } else if (length(put) == 3) {
      if (mean(put[c(1, 3)]) != put[2]) {
        stop("Value for heterozygotes must be equal to the mean of homozygotes")
      } else {
        x[x == put[1]] <- "AA"
        x[x == put[2]] <- "AB"
        x[x == put[3]] <- "BB"
      }
    }
  }
  # Names of all genotypes present in the marker data matrix.
  allgeno <- na.exclude(unique(c(x)))

  if (!all(nchar(allgeno) == 2)) {
    stop(cat("Alleles must be encoded by two letters."))
  }
  if (any(allgeno %in% c("XX", "XY", "YY"))) {
    stop("Letters 'X' and 'Y' are reserved for output")
  }
  if (mafThresh >= 1) {
    stop("A MAF threshold >= 1 will discard all loci.")
  }
  if (mafThresh < 0) {
    stop("The MAF threshold has to be positiive.")
  }

  output <- match.arg(output)
  # Separate heterozygous genotypes.
  namesplit <- strsplit(allgeno, split = "")
  heteroz <- allgeno[unlist(lapply(namesplit, function(x) x[1] != x[2]))]
  # Separate homozygous genotypes.
  homoz <- allgeno[!allgeno %in% c(missing, heteroz)]
  # Matrix for storage of genotype frequencies.
  genomatrix <- matrix(NA_real_, nrow = length(allgeno), ncol = ncol(x))
  rownames(genomatrix) <- allgeno
  colnames(genomatrix) <- colnames(x)
  # Frequency matrix of genotypes by marker.
  for (i in allgeno) {genomatrix[i, ] <- colSums(x == i)}
  # Genotypes excluding missings ("??").
  excna <- genomatrix[rownames(genomatrix) != missing, ]
  # Remove all markers with merely missing values.
  nona <- excna[, colSums(excna) != 0]
  # Remove all monomorphic markers.
  poly <- nona[, colSums(nona != 0) > 1]
  # Reduce the dimensions of the original data frame so that only polymorphic
  # markers are included.
  polyX <- x[, colnames(x) %in% colnames(poly)]
  # Reduce the genotype matrix to homozygous genotypes.
  homozmat <- poly[homoz, ]
  if (any(colSums(homozmat != 0) > 2)) stop(paste("markers with more than",
                                                  "two homozygous genotypes",
                                                  "present"))
  # For each locus, select the major and the minor genotype.
  n <- nrow(homozmat)
  major <- vapply(X = seq_len(ncol(homozmat)), FUN = function(i) {
    x <- homozmat[, i]
    majorGeno <- names(x[x == sort(x, partial = n)[n]])
    if(length(majorGeno) != 1) {
      majorGeno <- majorGeno[1]
    }
    return(majorGeno)
  }, FUN.VALUE = "character")

  minor <- vapply(X = seq_len(ncol(homozmat)), FUN = function(i) {
    x <- homozmat[,i]
    minorGeno <- names(x[x == sort(x, partial = n - 1)[n - 1]])
    if (length(minorGeno) != 1) {
      minorGeno <- minorGeno[2]
    }
    return(minorGeno)
  }, FUN.VALUE = "character")
  # For which loci are there no true major and minor genotypes?
  stopifnot(all(major != minor))

  majorMat <- matrix(major, nrow = 1)[rep(1, nrow(polyX)), ]
  minorMat <- matrix(minor, nrow = 1)[rep(1, nrow(polyX)), ]

  # Augment the matrix with reference/minor genotypes to have the same number
  # of rows as the input (Dent or Flint). Replace genotypes by new
  # identifiers for major, minor and heterozygous genotypes, respectively.
  polyX[polyX == majorMat] <- "XX"
  polyX[polyX == minorMat] <- "YY"
  polyX[!polyX %in% c("XX", "YY", missing)] <- "XY"
  n_maj <- colSums(polyX == "XX")
  n_min <- colSums(polyX == "YY")
  n_het <- colSums(polyX == "XY")
  # Number of genotypes, which are not missing.
  n_geno <- n_maj + n_min + n_het
  # Frequency of the major allele at any locus.
  p_maj <- ((2 * n_maj) + n_het) / (2 * n_geno)
  # Frequency of the minor allele at any locus.
  p_min <- 1 - p_maj
  stopifnot(max(p_min) <= 0.5)

  if (mafThresh == 0) {
    markernames <- colnames(polyX)
  } else {
    markernames <- colnames(polyX)[p_min >= mafThresh]
  }

  switch(output,
         markerNames = markernames,
         markerMAF = p_min,
         genoList = list(major_allele = major, minor_allele = minor))
}
