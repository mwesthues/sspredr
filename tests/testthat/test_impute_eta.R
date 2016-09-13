data("hybrid_nms", "mrna", "imp_snps")
geno <- vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 1,
               FUN.VALUE = character(1))
x <- imp_snps[rownames(imp_snps) %in% geno, ]
y <- mrna[rownames(mrna) %in% geno, ]
eta <- impute_eta(x = x, y = y, geno = geno, bglr_model = "BRR")

context("Content of ETA")
test_that("all random and fixed effects are included", {
  expect_length(eta, n = 3)
  expect_identical(vapply(eta, FUN = function(x) x$model,
                          FUN.VALUE = character(1)),
                   expected = c("FIXED", "BRR", "BRR"))
})

test_that("output dimensions match input", {
  expect_equal(ncol(eta[[1]][["X"]]), 2)
  expect_equal(ncol(eta[[2]][["X"]]), ncol(y))
})
