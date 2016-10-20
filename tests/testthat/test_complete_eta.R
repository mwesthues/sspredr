data("hybrid_nms", "imp_snps")
geno <- vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 1,
               FUN.VALUE = character(1))
x <- imp_snps[rownames(imp_snps) %in% geno, ]
eta_feature <- complete_eta(x = x, geno = geno, bglr_model = "BRR")
eta_kernel <- complete_eta(x = x, as_kernel = TRUE, geno = geno,
                           bglr_model = "BRR")

context("Content of ETA")
test_that("effects and model as expected", {
  expect_length(eta_feature, n = 1)
  expect_length(eta_kernel, n = 1)
  expect_identical(vapply(eta_feature, FUN = function(x) x$model,
                          FUN.VALUE = character(1)),
                   expected = "BRR")
  expect_identical(vapply(eta_kernel, FUN = function(x) x$model,
                          FUN.VALUE = character(1)),
                   expected = "BRR")
})

test_that("output dimensions match input", {
  expect_equal(ncol(eta_feature[[1]][["X"]]), ncol(x))
  expect_equal(ncol(eta_kernel[[1]][["X"]]), nrow(x))
  expect_equal(nrow(eta_feature[[1]][["X"]]), length(geno))
  expect_equal(nrow(eta_kernel[[1]][["X"]]), length(geno))
})
