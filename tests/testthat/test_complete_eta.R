data("hybrid_nms", "imp_snps")
geno <- vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 1,
               FUN.VALUE = character(1))
x <- imp_snps[rownames(imp_snps) %in% geno, ]
eta <- complete_eta(x = x, geno = geno, bglr_model = "BRR")

context("Content of ETA")
test_that("effects and model as expected", {
  expect_length(eta, n = 1)
  expect_identical(vapply(eta, FUN = function(x) x$model,
                          FUN.VALUE = character(1)),
                   expected = "BRR")
})

test_that("output dimensions match input", {
  expect_equal(ncol(eta[[1]][["X"]]), length(unique(geno)))
  expect_equal(nrow(eta[[1]][["X"]]), length(geno))
})
