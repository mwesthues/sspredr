context("Heterozygosity input")
test_that("no missing values allowed", {
  x <- imp_snps
  x_dim <- prod(dim(x))
  x_na <- sample(x_dim, size = 5, replace = FALSE)
  x[x_na] <- NA_real_
  expect_error(compute_het(x, output = "markerNames"),
               info = "NAs not allowed in 'x'")
})

test_that("heterozygosity exists", {
  x <- matrix(sample(c(0, 1), size = 200, replace = TRUE),
              nrow = 10, ncol = 20,
              dimnames = list(paste0("geno", seq_len(10)),
                              paste0("snp", seq_len(20))))
  expect_error(compute_het(x, output = "markerNames"))
})


context("Heterozygosity output")
test_that("output has correct class", {
  expect_is(compute_het(imp_snps, output = "markerNames"),
            class = "character")
  expect_is(compute_het(imp_snps, output = "markerHeterozygosity"),
            class = "numeric")
})


