# Generate transcriptomic data
data("hybrid_nms")
dent <- unique(vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 1,
                      FUN.VALUE = character(1)))
dent <- dent[seq_len(length(dent) * 0.5)]
flint <- unique(vapply(strsplit(hybrid_nms, split = "_"), FUN = "[[", 2,
                       FUN.VALUE = character(1)))
flint <- flint[seq_len(length(flint) * 0.5)]
set.seed(9034)
mrna <- matrix(rnorm(length(c(dent, flint)) * 500, mean = 0, sd = 5),
               nrow = length(c(dent, flint)), ncol = 500,
               dimnames = list(c(dent, flint),
                               paste0("col", seq_len(500))))
devtools::use_data(mrna, overwrite = TRUE)
