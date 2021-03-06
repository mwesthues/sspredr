% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/replace_na_with_major_genotype.R
\name{replace_na_with_major_genotype}
\alias{replace_na_with_major_genotype}
\title{Replace Missing Marker Genotypes}
\usage{
replace_na_with_major_genotype(dat, missing_value, major_genotype)
}
\arguments{
\item{dat}{A matrix. Genotype names are stored in rows whereas marker names
are stored in columns.}

\item{missing_value}{Specify the encoding of missing genotypes.}

\item{major_genotype}{Vector with the major allele at each locus.}
}
\value{
A matrix of the same type as \code{dat} but with missing values
 replaced by the major genotype at the corresponding locus.
}
\description{
\code{replace_na_with_major_genotype} replaces missing values with the most
 frequent genotype at each marker locus.
}
\examples{
 # Load a matrix with SNP genotypes encoded as numeric values
 data(marker_numeric)

 # Compute the major genotype at each locus.
 major_genotype <- sspredr::compute_maf(marker_numeric, output = "geno_list",
                                        missing = NA_real_,
                                        maf_threshold = 0)[["major_genotype"]]

 # Replace all missing genotypes with the major allele.
 replace_na_with_major_genotype(marker_numeric, missing_value = NA_real_,
                                major_genotype = major_genotype)
}
