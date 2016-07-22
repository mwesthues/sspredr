# Generate marker data with missing values.
set.seed(119)
x <- matrix(sample(c(0, 1, 2, NA_real_), size = 100,
                   prob = c(0.3, 0.1, 0.55, 0.05), replace = TRUE),
            nrow = 10, ncol = 10,
            dimnames = list(paste0("row", seq_len(10)),
                            paste0("col", seq_len(10))))
y <- x
storage.mode(y) <- "character"
y[y %in% "0"] <- "AA"
y[y %in% "1"] <- "AB"
y[y %in% "2"] <- "BB"
y[is.na(y)]  <- "??"
identical(is.na(x), y == "??")
marker_numeric <- x
marker_character <- y
devtools::use_data(marker_numeric, overwrite = TRUE)
devtools::use_data(marker_character, overwrite = TRUE)

# Generate marker data without missing values.
set.seed(9034)
imp_snps <- matrix(sample(c(0, 1, 2), size = 200,
                   prob = c(0.3, 0.05, 0.65), replace = TRUE),
                   nrow = 10, ncol = 20,
                   dimnames = list(paste0("row", seq_len(10)),
                                   paste0("col", seq_len(10))))
devtools::use_data(imp_snps, overwrite = TRUE)
