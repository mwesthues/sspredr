context("MAF input")
test_that("arguments play nice", {
  expect_error(compute_maf(marker_numeric, output = "marker_names",
                           missing = "??",
                           maf_threshold = 0),
               regexp = "Types of 'missing_value' and 'x' must be the same")
  expect_error(compute_maf(marker_character, output = "marker_names",
                           missing = NA_real_,
                           maf_threshold = 0),
             regexp = "Types of 'missing_value' and 'x' must be the same")
})


context("MAF recoding")
dummy_input <- matrix(
  c(
    2, 0, 0, 0, 2,
    0, NA_real_, 2, 2, 2,
    0, 0, 0, 2, 0
  ),
  ncol = 3
)

dummy_output <- matrix(
  c(
    0, 2, 2, 2, 0,
    0, NA_real_, 2, 2, 2,
    2, 2, 2, 0, 2
  ),
  ncol = 3
)

test_that("Major alles are defined correctly", {
  expect_equal(compute_maf(dummy_input, output = "recoded",
                           missing = NA_real_,
                           maf_threshold = 0),
               dummy_output)
})

context("MAF output")
test_that("output structure is correct", {
  expect_length(compute_maf(marker_numeric, output = "marker_names",
                            missing = NA_real_,
                            maf_threshold = 0),
                n = ncol(marker_numeric))
})

