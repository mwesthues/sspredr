# Using the mice data set from the 'BGLR' package, extract the numerator
# relationship matrix, a genomic relationship matrix and generate a
# transcriptomic relationship matrix from a multivariate normal distribution
# using the pedigree relationships to specify the covariance matrix of the
# mRNA vectors.
# Extract one phenotypic trait from the mice data set.
data("mice", package = "BGLR")
set.seed(943)
ped_nms <- sample(rownames(mice.A), size = 120)
mice_ped <- mice.A[match(ped_nms, rownames(mice.A)), 
                   match(ped_nms, colnames(mice.A))]
snp_nms <- sample(rownames(mice_ped), size = nrow(mice_ped) * 0.8,
                  replace = FALSE)
mice_snp <- mice.X[match(snp_nms, rownames(mice.X)), 
                   sample(colnames(mice.X), size = 600)]
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
