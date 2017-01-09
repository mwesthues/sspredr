data("mice_ped", "mice_snp", "mice_mrna")
library("magrittr")
library("purrr")
ped <- mice_ped
snp <- mice_snp
mrna <- mice_mrna
geno <- rownames(mice_ped)
bglr_model <- "BRR"
eta <- impute2(ped = mice_ped, snp = mice_snp, mrna = mice_mrna, geno = geno,
               bglr_model = bglr_model)
context("Content of ETA when imputing with pedigree and genomic info")
test_that("all random and fixed effects are included", {
  expect_length(eta, n = 4)
  expect_identical(vapply(eta, FUN = function(x) x$model,
                          FUN.VALUE = character(1)),
                   expected = c("FIXED", rep(get("bglr_model"), times = 3)))
})
