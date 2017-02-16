context("MAF input")
test_that("arguments play nice", {
  expect_error(compute_maf(marker_numeric, output = "marker_names",
                           missing = "??",
                           maf_threshold = 0),
               regexp = "Types of 'missing_value' and 'x' must be the same")
})


context("MAF output")
test_that("output structure is correct", {
  expect_length(compute_maf(marker_numeric, output = "marker_names",
                            missing = NA_real_,
                            maf_threshold = 0),
                n = ncol(marker_numeric))
})

test_that("Length of major and minor genotypes are the same", {
  multiple_genotypes <- matrix(c(
   "AA", "AA", "AT", "TT",
   "CC", "CG", "GG", "GG",
   "TT", "NN", "TT", "AA",
   "CC", "CC", "CC", "CC"
  ), ncol = 4, dimnames = list(NULL, paste0("col", seq_len(4))))
  maf_geno_length <- compute_maf(multiple_genotypes, output = "geno_list",
                                 missing = "NN",
                                 maf_threshold = 0)
  expect_identical(length(maf_geno_length$major_genotype),
                   length(maf_geno_length$minor_genotype))
})

test_that("Major and minor genotypes are unique", {
  geno_list <- compute_maf(marker_character,
                           output = "geno_list",
                           missing = "??",
                           maf_threshold = 0.9)
  expect_false(any(geno_list$major_genotype == geno_list$minor_genotype))
})


test_that("Input without heterozygous genotypes is accepted", {
  homozygous <- matrix(sample(c(0, 1), size = 100, replace = TRUE),
                       nrow = 10, ncol = 10)
  colnames(homozygous) <- paste0("col", seq_len(10))
  expect_silent(compute_maf(homozygous,
                            output = "marker_names",
                            missing = NA_real_,
                            maf_threshold = 0))
})
