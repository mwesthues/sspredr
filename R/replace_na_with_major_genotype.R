#' Replace Missing Marker Genotypes
#'
#' \code{replace_na_with_major_genotype} replaces missing values with the most
#'  frequent genotype at each marker locus.
#'
#' @param dat A matrix. Genotype names are stored in rows whereas marker names
#'  are stored in columns.
#' @param missing_value Specify the encoding of missing genotypes.
#' @param major_genotype Vector with the major allele at each locus.
#' @return A matrix of the same type as \code{dat} but with missing values
#'  replaced by the major genotype at the corresponding locus.
#' @examples
#'  # Load a matrix with SNP genotypes encoded as numeric values
#'  data(marker_numeric)
#'
#'  # Compute the major genotype at each locus.
#'  major_genotype <- sspredr::compute_maf(marker_numeric, output = "geno_list",
#'                                         missing = NA_real_,
#'                                         maf_threshold = 0)[["major"]]
#'
#'  # Replace all missing genotypes with the major allele.
#'  replace_na_with_major_genotype(marker_numeric, missing_value = NA_real_,
#'                                 major_genotype = major_genotype)
#' @importFrom magrittr %>%
#' @export
replace_na_with_major_genotype <- function(dat, missing_value, major_genotype) {

  if (typeof(dat) != typeof(missing_value)) {
    stop("'missing_value' must have the same type as the entries in 'dat'")
  }

  if (typeof(dat) %in% c("integer", "double") &&
      !typeof(major_genotype) %in% c("integer", "double")) {
    stop("'major_genotype' must have the same type as the entries in 'dat'")
  }

  no_na_snp <- dat %>%
    ncol() %>%
    seq_len() %>%
    lapply(function(i) {
      x <- dat[, i]
      if (anyNA(x)) {
        x[is.na(x)] <- major_genotype[i]
      } else if (!is.na(missing_value) && any(x == missing_value)) {
        x[x == missing_value] <- major_genotype[i]
      }
      unique_genotypes <- x %>%
        unique() %>%
        length()
      x
    }) %>%
    unlist() %>%
    matrix(., ncol = ncol(dat), byrow = FALSE)
  rownames(no_na_snp) <- rownames(dat)
  colnames(no_na_snp) <- colnames(dat)
  no_na_snp
}
