data("hybrid_nms", "mrna", "imp_snps")
geno <- vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 1,
               FUN.VALUE = character(1))
x <- imp_snps[rownames(imp_snps) %in% geno, ]
y <- mrna[rownames(mrna) %in% geno, ]
eta_kernel <- impute_eta(x = x, y = y, geno = geno, as_kernel = TRUE,
                         is_pedigree = FALSE, bglr_model = "BRR")

context("Content of ETA when imputing only with genomic info")
test_that("all random and fixed effects are included", {
  expect_length(eta_kernel, n = 3)
  expect_identical(vapply(eta_kernel, FUN = function(x) x$model,
                          FUN.VALUE = character(1)),
                   expected = c("FIXED", "BRR", "BRR"))
})

