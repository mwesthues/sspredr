# Using the mice data set from the 'BGLR' package, extract the numerator
# relationship matrix, a genomic relationship matrix and generate a
# transcriptomic relationship matrix from a multivariate normal distribution
# using the pedigree relationships to specify the covariance matrix of the
# mRNA vectors.
# Extract one phenotypic trait from the mice data set.
data("mice", package = "BGLR")
set.seed(30943)
mice_ped <- mice.A[seq_len(150), seq_len(150)]
mice_snp <- build_kernel(mice.X[seq_len(150), ])
snp_nms <- sample(rownames(mice_ped), size = nrow(mice_ped) * 0.8,
                  replace = FALSE)
mice_snp <- mice_snp[match(snp_nms, rownames(mice_snp)),
                     match(snp_nms, colnames(mice_snp))]
mrna_nms <- sample(snp_nms, size = nrow(mice_snp) * 0.8, replace = FALSE)
mu <- rnorm(n = nrow(mice_ped), mean = 10, sd = 2)
mice_mrna <- t(MASS::mvrnorm(n = 400, mu = mu, Sigma = mice_ped))
mice_mrna <- mice_mrna[match(mrna_nms, rownames(mice_mrna)), ]
mice_pheno <- mice.pheno[mice.pheno$SUBJECT.NAME %in% rownames(mice_ped),
                         c("SUBJECT.NAME", "Obesity.BMI")]
devtools::use_data(mice_ped, overwrite = TRUE)
devtools::use_data(mice_snp, overwrite = TRUE)
devtools::use_data(mice_mrna, overwrite = TRUE)
devtools::use_data(mice_pheno, overwrite = TRUE)
