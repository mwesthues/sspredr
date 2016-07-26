#' SNP genotypes (character)
#'
#' A dataset containing the genotypes of ten individuals for ten loci. The
#' variables are as follows:
#'
#' @format A matrix with 10 rows and 10 columns. All genotypes have class
#' 'character'. The matrix contains missing values.
"marker_character"


#' SNP genotypes (numeric)
#'
#' A dataset containing the genotypes of ten individuals for ten loci. The
#' variables are as follows:
#'
#' @format A matrix with 10 rows and 10 columns. All genotypes have class
#' 'numeric'. The matrix contains missing values.
"marker_numeric"


#' Imputed, numeric SNP genotypes
#'
#' A dataset containing the genotypes of ten individuals for 20 loci. The
#' variables are as follows:
#'
#' @format A matrix with 10 rows and 20 columns. All genotypes have class
#' 'numeric'. No missing values exist in the data set.
"imp_snps"

#' Hybrid names.
#'
#' A vector containing the names of hybrids. It can be used for sampling a
#' cross-validation scheme.
#'
#' @format A character vector with hybrid names where the two parents are
#'  separated by an underscore '_'.
"hybrid_nms"
