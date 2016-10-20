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


#' Transcriptomic data
#'
#' A matrix containing the gene expression values for 44 genotypes and 500
#' mRNAs.
#'
#' @format A matrix with 44 rows and 500 columns.
"mrna"

#' Phenotypic data
#'
#' A matrix containing the agronomic data for 532 hybrids and two traits.
#'
#' @format A matrix with 532 rows and two columns.
"Pheno"


#' Mice pedigree data
#'
#' A numerator relationship matrix with pedigree records for 150 mice.
#'
#' @format A matrix with 150 rows and 150 columns.
"mice_ped"


#' Mice genomic relationship data
#'
#' A genomic relationship matrix with records for 120 mice.
#'
#' @format A matrix with 120 rows and 120 columns.
"mice_snp"


#' Mice transcriptomic data
#'
#' A matrix with transcriptomic features for 96 mice. The values were calculated
#' using a multivariate normal distribution where pedigree records in 'mice_ped'
#' specified the covariance matrix of the mRNA vectors.
#'
#' @format A matrix with 96 rows and 400 columns.
"mice_mrna"


#' Mice phenotypic data
#'
#' A data frame with obesity records for 150 mice taken from the 'BGLR' package.
#'
#' @format A data.frame with 150 rows and 2 columns.
"mice_pheno"

