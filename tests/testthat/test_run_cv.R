data("hybrid_nms", "mrna", "imp_snps", "Pheno")
geno <- vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 1,
               FUN.VALUE = character(1))
x <- imp_snps[rownames(imp_snps) %in% geno, ]
y <- mrna[rownames(mrna) %in% geno, ]
# BGLR ETA objects
eta <- impute_eta(x = x, y = y, geno = geno, as_kernel = TRUE,
                  is_pedigree = FALSE, bglr_model = "BRR")
eta[] <- lapply(seq_along(eta), FUN = function(i) {
  dat <- eta[[i]]
  x <- dat[["X"]]
  rownames(x) <- hybrid_nms
  dat[["X"]] <- x
  dat
})
# CV scheme generation
cv_mat <- sample_cv(hybrid_nms, n_mother = 39, n_father = 33,
                    n_hyb_trn = 200, min_size = 20, rounds = 2L,
                    hybrid_split = "_")
cv_mat <- cv_mat[[1]]
colnames(cv_mat) <- seq_len(ncol(cv_mat))
rownames(cv_mat) <- gsub("DF_", replacement = "", x = rownames(cv_mat))
cv <- as.data.frame(cv_mat)
cv$Sample_ID <- rownames(cv)
cv <- tidyr::gather(cv, Run, Set, -Sample_ID)
dir.create("./bglr_out")
trait <- "gtm"
run <- 1L
father_idx <- 2
mother_idx <- 1
split_char <- "_"
father <- unlist(strsplit(hybrid_nms[run], split = split_char))[father_idx]
mother <- unlist(strsplit(hybrid_nms[run], split = split_char))[mother_idx]
pred_res <- run_cv(Pheno = Pheno,
                   ETA = eta,
                   cv = cv,
                   father_idx = father_idx,
                   mother_idx = mother_idx,
                   split_char = split_char,
                   trait = trait,
                   iter = 2500L,
                   speed_tst = FALSE,
                   run = run,
                   verbose = FALSE,
                   out_loc = "./bglr_out/")
speed_res <- run_cv(Pheno = Pheno,
                    ETA = eta,
                    cv = cv,
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
  expect_equal(run, unique(pred_res$Run))
})

context("Speed test")
test_that("speed test worked", {
  expect_is(speed_res, class = "proc_time")
})
