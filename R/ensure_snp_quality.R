#' Ensure SNP Quality
#'
#' \code{ensure_snp_quality} applies quality checks to SNP data
#'
#' @param snp A matrix. Genotype names are stored in rows whereas marker names
#'  are stored in columns.
#' @param callfreq_check Logical scalar. Shall markers with a low call
#'  frequency (*i.e.* a high number of missing marker genotypes) be removed?
#' @param callfreq_threshold Numeric scalar. Minimum level of non-missing
#'  genotypes per locus.
#' @param maf_check Logical scalar. Shall markers with a low minor allele
#'  frequency be removed?
#' @param maf_threshold Numeric scalar. Minimum frequency of the minor allele
#'  at each locus.
#' @param any_missing Logical scalar. Does the input matrix \code{snp} contain
#'  any missing values? If this is the case, they will be replaced with the
#'  major genotype at this locus.
#' @param missing_value Specify the encoding of missing genotypes.
#' @param remove_duplicated Logical scalar. Should only unique marker loci be
#'  returned?
#' @return Depending on the choice of parameters \code{callfreq_check} and
#'  \code{maf_check}, respectively, \code{ensure_snp_quality} will return a
#'  matrix with marker genotypes that have passed important quality checks.
#' @examples
#'  # Load a matrix with SNP genotypes encoded as numeric values
#'  data(marker_numeric)
#'
#'  # Add a duplicated marker locus to the data.
#'  snp <- marker_numeric
#'  snp <- cbind(snp, snp[, 1])
#'
#'  # Return markers without missing values, a call frequency greater or equal
#'  # to 0.95 and a minor allele frequency greater or equal to 0.6. Finally,
#'  # remove all duplicated markers.
#'  ensure_snp_quality(snp, maf_threshold = 0.1, missing_value = NA_real_)
#'
#'  # Load a matrix with SNP genotypes encoded as character values
#'  data(marker_character)
#'
#'  # Return markers without missing values and a call frequency equal to or
#'  # greater than 0.9. Keep duplicated markers if present.
#'  ensure_snp_quality(marker_character, callfreq_threshold = 0.9,
#'                     maf_check = FALSE, missing_value = "??",
#'                     remove_duplicated = FALSE)
#' @importFrom magrittr %>%
#' @export
ensure_snp_quality <- function(snp, callfreq_check = TRUE,
                               callfreq_threshold = 0.95, maf_check = TRUE,
                               maf_threshold = 0.05, any_missing = TRUE,
                               missing_value = NA_character_,
                               remove_duplicated = TRUE) {

  if (!class(snp) == "matrix") stop("'snp' is not of class 'matrix'")
  if (isTRUE(callfreq_check)) {
    if (!is.numeric(callfreq_threshold)) {
      stop("Specify a call frequency threshold")
    }
    if (callfreq_threshold < 0 || callfreq_threshold > 1) {
      stop("Call frequency must lie in within the range from 0 to 1")
    }

    passed_cf_nms <- sspredr::compute_cf(snp, output = "marker_names",
                                         missing = missing_value,
                                         call_threshold = callfreq_threshold)
    snp <- snp[, colnames(snp) %in% passed_cf_nms]
  }

  if (isTRUE(maf_check)) {
    if (!is.numeric(maf_threshold)) {
      stop("Specify a minor allele frequency threshold")
    }
    if (maf_threshold < 0 || maf_threshold > 1) {
      stop("Minor allele frequency must lie in within the range from 0 to 1")
    }

    passed_maf_nms <- sspredr::compute_maf(snp,
                                           output = "marker_names",
                                           missing = missing_value,
                                           maf_threshold = maf_threshold)
    snp <- snp[, colnames(snp) %in% passed_maf_nms]
  }

  if (isTRUE(any_missing)) {
    major_genotype <- snp %>%
      sspredr::compute_maf(
      output = "geno_list", missing = missing_value, maf_threshold = 0
      ) %>%
      .[["major"]]
      snp <- replace_na_with_major_genotype(snp, missing_value = missing_value,
                                            major_genotype = major_genotype)
  }

  if (isTRUE(remove_duplicated)) {
    snp <- unique(snp, MARGIN = 2)
  }
  snp
}
