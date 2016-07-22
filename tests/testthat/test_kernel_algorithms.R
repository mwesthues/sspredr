data(mice.X)
mice.snp <- mice.X[1:200, ]
char_snps <- imp_snps
storage.mode(char_snps) <- "character"
char_snps[char_snps == "0"] <- "BB"
char_snps[char_snps == "1"] <- "AB"
char_snps[char_snps == "2"] <- "AA"

test_that("output has dimension equal to number of individuals", {
  expect_identical(prod(dim(build_kernel(mice.snp))), nrow(mice.snp)^2)
})

test_that("resulting relationship matrix can be inverted successfully", {
  expect_error(solve(build_kernel(mice.snp, algorithm = "RadenII")), NA)
  expect_error(solve(build_kernel(mice.snp, algorithm = "Zhang")), NA)
})

test_that("only numeric input is accepted", {
  expect_error(build_kernel(char_snps))
  expect_error(build_kernel(imp_snps), NA)
})
