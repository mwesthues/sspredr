char_geno_lst <- compute_maf(marker_character, output = "genoList")
char_major <- char_geno_lst[["major_allele"]]
char_minor <- char_geno_lst[["minor_allele"]]

num_geno_lst <- compute_maf(marker_numeric, output = "genoList",
                            missing = NA_real_)
num_major <- num_geno_lst[["major_allele"]]
num_minor <- num_geno_lst[["minor_allele"]]



context("Encoding")
test_that("encoding is correct", {
  expect_error(recode_snps(marker_character,
                           major = char_major,
                           minor = char_minor,
                           major_coding = 1,
                           minor_coding = 0,
                           het_coding = 2,
                           na_coding = NA_real_))
  expect_error(recode_snps(marker_numeric,
                           major = num_major,
                           minor = num_minor,
                           major_coding = "AA",
                           minor_coding = "BB",
                           het_coding = "CC",
                           na_coding = NA_character_))
  expect_error(recode_snps(marker_character, major = char_major,
                           minor = char_minor,
                         major_coding = 1, minor_coding = "0", het_coding = 2,
                         na_coding = NA_real_),
             info = "classes of encodings are not identical")
})

test_that("class of ouptput equal to class of encoding", {
  expect_is(c(recode_snps(marker_character,
                          major = char_major,
                          minor = char_minor,
                          major_coding = 1,
                          minor_coding = 0,
                          het_coding = 0.5,
                          na_coding = NA_real_)),
            class = "numeric")
  expect_is(c(recode_snps(marker_numeric,
                          major = num_major,
                          minor = num_minor,
                          major_coding = "AA",
                          minor_coding = "BB",
                          het_coding = "AB",
                          na_coding = NA_character_)),
            class = "character")
})

context("Recode output")
test_that("more SNPs than genotypes used", {
   expect_error(recode_snps(marker_character[, 1:6],
                            major = char_major,
                            minor = char_minor,
                            major_coding = 2, minor_coding = 0, het_coding = 1,
                            na_coding = NA_real_),
               regexp = "non-conformable arrays",
               info = "Less SNP-columns than genotype-rows")
})

test_that("output has correct class", {
  expect_is(recode_snps(marker_character,
                       major = char_major,
                       minor = char_minor,
                       major_coding = 1,
                       minor_coding = 0,
                       het_coding = 0.5,
                       na_coding = NA_real_),
            class = "matrix")
})

test_that("output has same dimensions as input", {
  expect_identical(dim(recode_snps(marker_character,
                                   major = char_major,
                                   minor = char_minor,
                                   major_coding = 1,
                                   minor_coding = 0,
                                   het_coding = 0.5,
                                   na_coding = NA_real_)),
                   expected = dim(marker_character))
})


test_that("conversion of genotypes succeeded", {
  expect_true(all(c(recode_snps(marker_numeric,
                                major = num_major,
                                minor = num_minor,
                                major_coding = "AA",
                                minor_coding = "BB",
                                het_coding = "AB",
                                na_coding = NA_character_)) %in%
                    c("AA", "BB", "AB", NA_character_)))
})





