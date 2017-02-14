#' Minor Allele Frequency
#'
#' \code{compute_maf} returns results from minor allele frequency checks
#'
#' @param x A matrix. Genotype names are stored in rows whereas marker names
#'  are stored in columns.
#' @param output Character vector with output options.
#' @param missing_value Character vector providing the encoding of missing elements.
#' @param maf_threshold Numeric or complex vector specifying the minor allele
#'  frequency threshold.
#' @return If \code{output} is "marker_names" a character vector with marker
#'  names that have passed the quality check will be returned.
#'  If \code{output} is "marker_maf" a numeric vector with the minor allele
#'  frequency at each marker locus will be returned.
#'  If \code{output} is "geno_list" a list with encoding for the major and
#'  the minor genotype at each locus will be returned.
#' @examples
#'  # Load a matrix with SNP genotypes encoded as numeric values
#'  data(marker_numeric)
#'
#'  # Return the names of all polymorphic markers.
#'  compute_maf(marker_numeric, output = "marker_names",
#'              missing_value = NA_real_, maf_threshold = 0)
#'
#'  # Set one locus to monomorphic
#'  marker_mono <- marker_numeric
#'  marker_mono[, 1] <- 1
#'  compute_maf(marker_numeric, output = "marker_names",
#'              missing_value = NA_real_, maf_threshold = 0)
#'
#'  # Load a matrix with SNP genotypes encoded as character values
#'  data(marker_character)
#'
#'  # Return the minor allele frequency at each locus with MAF >= 0.05.
#'  compute_maf(marker_character, output = "marker_maf",
#'              missing_value = "??",
#'              maf_threshold = 0.05)
#'
#'  # Return the minor and major genotype at each locus.
#'  compute_maf(marker_character, output = "geno_list", missing_value = "??",
#'              maf_threshold = 0)
#'
#' # Examine a case where there more than three genotypes.
#' multiple_genotypes <- matrix(c(
#'  "AA", "AA", "AT", "TT",
#'  "CC", "CG", "GG", "GG",
#'  "TT", "NN", "TT", "AA",
#'  "CC", "CC", "CC", "CC"
#' ), ncol = 4, dimnames = list(NULL, paste0("col", seq_len(4))))
#' compute_maf(multiple_genotypes, output = "marker_names",
#'             missing_value = "NN", maf_threshold = 0)
#' @export
compute_maf <- function(
  x,
  output = c("marker_names", "marker_maf", "recoded", "geno_list"),
  missing_value = NA_character_,
  maf_threshold = 0) {

  # Store the original type of the input so that results for the
  # 'output'-option 'recoded' can be converted back to the original type.
  x_type <- typeof(x)

  # General input checks
  if (class(x) != "matrix") {
    stop("Input has to be of class 'matrix'")
  }

  if (typeof(missing_value) != typeof(x) && !is.na(missing_value)) {
    stop("Types of 'missing_value' and 'x' must be the same.")
  }

  # Recode missing values as 'NA' so that they can easily be excluded by using
  # R-specific functions.
  if (!is.na(missing_value)) {
    x[x == missing_value] <- NA
  }

  # Extract all unique genotpyes that are present in the data. This is necessary
  # for computationally efficient counting of elements at each locus.
  unique_elements <- sort(stats::na.exclude(unique(c(x))))

  # Count the number of occurences of each possible genotype at any locus. Store
  # the results in a matrix with number of rows equal to the number of unique
  # genotypes and number of columns equal to the number of loci.
  count_lst <- lapply(unique_elements, FUN = function(y) {
    colSums(x == y, na.rm = TRUE)
  })
  names(count_lst) <- unique_elements
  count_mat <- do.call(rbind, count_lst)
  colnames(count_mat) <- colnames(x)

  # Remove all monomorphic loci. First, they are uninformative and will only
  # increase noise in the models. Second, this practice makes it easier to
  # determine the minor genotype at each locus by looking for the smallest value
  # that is non-zero. Otherwise, the algorithm would always pick a genotype from
  # 'count_mat' that does not even exist at this particular locus.
  polymorphic_loci <- vapply(seq_len(ncol(count_mat)), FUN = function(i) {
    non_zero_count <- sum(count_mat[, i] != 0)
    non_zero_count > 1 && non_zero_count != 0
    sum(count_mat[, i] != 0) > 1 && sum(count_mat[, i])
  }, FUN.VALUE = logical(1))
  count_mat <- count_mat[, polymorphic_loci]


  # Determine the hom genotypes at each locus. Likewise, determine the
  # het genotypes at each locus. Without this separation there would
  # likely be cases where the frequency of het genotypes is higher than
  # the frequency of the minor genotype. This would exacerbate an efficient
  # determination of the minor genotype.
  if (isTRUE(x_type == "character")) {
    hom_alleles <- vapply(strsplit(rownames(count_mat), split = ""),
                                 FUN = function(x) {
      x[1] == x[2]
    }, FUN.VALUE = logical(1))
    hom_count_mat <- count_mat[hom_alleles, , drop = FALSE]
    het_count_mat <- count_mat[!hom_alleles, , drop = FALSE]
  }

  if (length(unique_elements) == 2) {
    hom_count_mat <- count_mat
  } else if (length(unique_elements) == 3) {
    el1 <- unique_elements[1]
    el2 <- unique_elements[2]
    el3 <- unique_elements[3]
    if (x_type %in% c("double", "integer")) {
      hom_mean <- mean(c(el1, el3))
      stopifnot(hom_mean == el2)
    }
    het_count_mat <- count_mat[rownames(count_mat) == el2, , drop = FALSE]
    hom_count_mat <- count_mat[rownames(count_mat) != el2, , drop = FALSE]
  }
  # With the previously assembled objects for the determination of the type of
  # the genotype (major, minor, het), determine their frequencies and
  # also return the matches.
  allele_type_lst <- lapply(seq_len(ncol(count_mat)), FUN = function(i) {
    locus <- stats::na.exclude(x[, i])
    major <- names(which.max(hom_count_mat[, i]))
    non_zero_hom <- hom_count_mat[, i] != 0
    minor <- names(which.min(hom_count_mat[non_zero_hom, i]))
    putative_het <- rownames(het_count_mat)
    het <- putative_het[het_count_mat[, i] != 0]
    n_het <- sum(locus == het)
    n_major <- sum(locus == major)
    n_minor <- sum(locus == minor)
    # Number of genotypes, which are not missing.
    n_geno <- n_major + n_het + n_minor
    # Frequency of the major allele at any locus.
    p_major <- ((2 * n_major) + n_het) / (2 * n_geno)
    # Frequency of the minor allele at any locus.
    p_minor <- 1 - p_major

    if (isTRUE(length(het) == 0)) het <- NA

    list(major_genotype = major,
         het_genotype = het,
         minor_genotype = minor,
         major_frequency = p_major,
         minor_frequency = p_minor)
  })
  names(allele_type_lst) <- colnames(count_mat)

  # Transform the results so that each component (major_genotype, het_genotype,
  # minor_genotype, major_frequency, minor_frequency) is returned as a separate
  # vector that can be easily extracted subsequently.
  transposed_lst <- purrr::transpose(allele_type_lst)
  res_lst <- lapply(transposed_lst, FUN = unlist)

  # Names of marker loci, which passed the minor allele frequency threshold.
  marker_names <- res_lst %>%
    .["minor_frequency"] %>%
    purrr::flatten_dbl() %>%
    purrr::keep(. >= maf_threshold) %>%
    names()

  # Return the names of the major and minor genotypes at each locus.
  geno_list <- res_lst %>%
    purrr::keep(names(.) %in% c("major_genotype", "minor_genotype"))
  geno_list[] <- lapply(geno_list, FUN = function(x) {
    storage.mode(x) <- x_type
    x
  })

  marker_maf <- res_lst %>%
    purrr::keep(names(.) %in% c("major_frequency", "minor_frequency"))

  output <- match.arg(output)
  switch(output,
         marker_names = marker_names,
         marker_maf = marker_maf,
         geno_list = geno_list)
}
