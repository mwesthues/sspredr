library("purrr")
library("magrittr")
data("marker_numeric")

context("Ensure SNP Quality")
# Add a marker locus with a minor allele frequency of 0.1.
dat <- cbind(marker_numeric, matrix(c(rep(2, times = 9), 1), ncol = 1))
# Add a duplicated locus.
dat <- cbind(dat, dat[, 5])
colnames(dat) <- paste0("col", seq_len(ncol(dat)))

test_that("Columns 3, 4, 11 and 12 get killed", {
  out1 <- sspredr::ensure_snp_quality(
    snp = dat,
    callfreq_check = TRUE,
    callfreq_threshold = 0.9,
    maf_check = TRUE,
    maf_threshold = 0.1,
    any_missing = TRUE,
    missing_value = NA_real_,
    remove_duplicated = TRUE
  )
  expect_equal(c("col1", "col2", "col5", "col6", "col7", "col8", "col9",
                 "col10"),
               colnames(out1))
})

test_that("Columns 4, 10 and 11 get killed", {
  out2 <- sspredr::ensure_snp_quality(
    snp = dat,
    callfreq_check = TRUE,
    callfreq_threshold = 0.8,
    maf_check = TRUE,
    maf_threshold = 0.3,
    any_missing = TRUE,
    missing_value = NA_real_,
    remove_duplicated = FALSE
  )
  expect_equal(c("col1", "col2", "col3", "col5", "col6", "col7", "col8", "col9",
                 "col12"),
               colnames(out2))
})
