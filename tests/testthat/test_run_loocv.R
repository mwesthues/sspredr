data("hybrid_nms", "mrna", "imp_snps", "Pheno")
geno <- vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 1,
               FUN.VALUE = character(1))
x <- imp_snps[rownames(imp_snps) %in% geno, ]
y <- mrna[rownames(mrna) %in% geno, ]
eta <- impute_eta(x = x, y = y, geno = geno, 
                  as_kernel = TRUE, is_pedigree = FALSE, bglr_model = "BRR")
eta[] <- lapply(seq_along(eta), FUN = function(i) {
  dat <- eta[[i]]
  x <- dat[["X"]]
  rownames(x) <- hybrid_nms
  dat[["X"]] <- x
  dat
})
dir.create("./bglr_out")
trait <- "gtm"
run <- 1L
father_idx <- 2
mother_idx <- 1
split_char <- "_"
father <- unlist(strsplit(hybrid_nms[run], split = split_char))[father_idx]
mother <- unlist(strsplit(hybrid_nms[run], split = split_char))[mother_idx]
pred_res <- run_loocv(Pheno = Pheno,
                      ETA = eta,
                      hybrid = TRUE,
                      father_idx = father_idx,
                      mother_idx = mother_idx,
                      split_char = split_char,
                      trait = trait,
                      iter = 2500L,
                      speed_tst = FALSE,
                      run = run,
                      verbose = FALSE,
                      out_loc = "./bglr_out/")
speed_res <- run_loocv(Pheno = Pheno,
                       ETA = eta,
                       hybrid = TRUE,
                       father_idx = father_idx,
                       mother_idx = mother_idx,
                       split_char = split_char,
                       trait = trait,
                       iter = 250L,
                       speed_tst = TRUE,
                       run = run,
                       verbose = FALSE,
                       out_loc = "./bglr_out/")
unlink("./bglr_out", recursive = TRUE)

context("Prediction")
test_that("prediction worked", {
  expect_is(pred_res, "data.frame")
  expect_false(anyNA(pred_res))
  expect_equal(pred_res$Geno, hybrid_nms[run])
  expect_equal(pred_res$Phenotype, trait)
})

test_that("hybrid prediction worked", {
  expect_equal(pred_res$Father, father)
  expect_equal(pred_res$Mother, mother)
})

context("Speed test")
test_that("speed test worked", {
  expect_is(speed_res, class = "proc_time")
})

