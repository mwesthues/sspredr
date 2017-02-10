library("purrr")
library("magrittr")
data("marker_numeric")

context("Replacement of missing values with major genotype")
geno_list <- sspredr::compute_maf(marker_numeric, output = "geno_list",
                                  missing_value = NA_real_, maf_threshold = 0)
major_genotypes <- geno_list$major_genotype
n_geno <- nrow(marker_numeric)
expected_major_geno_count <- vapply(seq_len(ncol(marker_numeric)),
                                    FUN = function(i) {
  locus <- marker_numeric[, i]
  major_count <- sum(locus == major_genotypes[i], na.rm = TRUE)
  na_count <- sum(is.na(locus))
  exp_future_major_count <- major_count + na_count
  exp_future_major_count
}, FUN.VALUE = integer(1))
imputed_snps <- sspredr::replace_na_with_major_genotype(
  marker_numeric,
  missing_value = NA_real_,
  major_genotype = major_genotypes
  )
realized_major_geno_count <- vapply(seq_len(ncol(imputed_snps)),
                                    FUN = function(i) {
  locus <- imputed_snps[, i]
  major_count <- sum(locus == major_genotypes[i])
  major_count
}, FUN.VALUE = integer(1))

test_that("missing values were replaced with the major genotype", {
  expect_equal(expected_major_geno_count, realized_major_geno_count)
})
