data("hybrid_nms")
# Generate phenotypic data.
set.seed(119)
n <- length(hybrid_nms)
gtm <- sample(rnorm(n, mean = 190, sd = 4), replace = FALSE)
gts <- sample(rnorm(n, mean = 0.25, sd = 0.03), replace = FALSE)
Pheno <- cbind(gtm, gts)
rownames(Pheno) <- hybrid_nms
devtools::use_data(Pheno, overwrite = TRUE)
