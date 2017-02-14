#' Recode genotypic data
#'
#' \code{recode_snps} recodes the elements of a SNP-matrix
#'
#' @param x A matrix. Genotype names are stored in rows whereas marker names
#' are stored in columns.
#' @param major Character vector with names of major genotype at each locus.
#' @param minor Character vector with names of minor genotype at each locus.
#' @param major_coding Scalar specifying how to recode the major genotype.
#' @param minor_coding Scalar specifying how to recode the minor genotype.
#' @param het_coding Scalar specifying how to recode heterozygous loci.
#' @param missing_value Scalar specifying how missing values are encoded in \code{x}.
#' @param na_coding Scalar specifying how to recode missing values.
#' @return A matrix with the same dimensions as \code{x} but with recoded SNP
#'  genotypes.
#' @examples
#'  # First, determine the major and minor allele at each locus of SNP-genotypes
#'  geno_lst <- compute_maf(marker_character, output = "geno_list",
#'                          missing_value = "??", maf_threshold = 0)
#'  major <- geno_lst[["major_genotype"]]
#'  minor <- geno_lst[["minor_genotype"]]
#'
#'  # Recode the major allele as '2', the minor allele as '0' and heterozygotes
#'  # as '0'.
#'  recode_snps(marker_character, major = major, minor = minor,
#'              major_coding = 2, minor_coding = 0, het_coding = 1,
#'              missing_value = NA_character, na_coding = NA_real_)
#' @export
recode_snps <- function(x, major, minor, major_coding, minor_coding,
                        het_coding, missing_value, na_coding) {
  if (class(x) != "matrix") {
    stop("'x' is not a matrix")
  }
  class_encoding <- unique(vapply(list(major_coding, minor_coding, het_coding,
                                       na_coding),
                                  FUN = class, FUN.VALUE = character(1)))
  if (length(class_encoding) != 1) {
    stop("the classes of all encodings are not identical")
  }
  if (class_encoding %in% c("integer", "numeric")) {
    if (mean(c(minor_coding, major_coding)) != het_coding) {
      stop("Value for heterozygotes not equal to the mean of homozygotes")
    }
  } else if (class_encoding == "character") {
    hom_encoding <- all(unlist(lapply(list(major_coding, minor_coding),
                                      FUN = function(x) {
      length(unique(unlist(strsplit(x, split = "")))) == 1
    })))
    if (!isTRUE(hom_encoding)) {
      stop("Homozygous genotypes encoded as heterozygous genotypes")
    }
    het_major <- grepl(stringr::str_sub(het_coding, start = 1, end = 1),
                       x = major_coding)
    het_minor <- grepl(stringr::str_sub(het_coding, start = 2, end = 2),
                       x = minor_coding)
    if (!all(c(het_major, het_minor))) {
      stop("'het_coding' does not match 'major_coding' and/or 'minor_coding'")
    }
  }


  major_mat <- matrix(major, nrow = 1)
  minor_mat <- matrix(minor, nrow = 1)
  x[x == major_mat[rep(1, nrow(x)), ]] <- major_coding
  x[x == minor_mat[rep(1, nrow(x)), ]] <- minor_coding
  if (is.na(missing_value)) {
    x[is.na(x)] <- na_coding
  } else {
    x[x == missing_value] <- na_coding
  }
  x[!x %in% c(get("major_coding"),
              get("minor_coding"), get("na_coding"))] <- het_coding
  storage.mode(x) <- class_encoding
  x
}
