test_that("'missing' has the same type as the elements of 'x'", {
  expect_error(compute_cf(marker_numeric, output = "markerNames",
                          missing = "??"))

})

test_that("'callThresh' is required when 'markerNames' are expected", {
  expect_error(compute_cf(marker_character, output = "markerNames",
                          missing = "??"))
})

test_that("output is correct", {
  expect_is(compute_cf(marker_character, output = "markerNames",
                       missing = "??", callThresh = 0.9),
            class = "character")
  expect_is(compute_cf(marker_numeric, output = "markerNames",
                       missing = NA_real_, callThresh = 0.9),
            class = "character")
  expect_length(compute_cf(marker_character, output = "markerNames",
                       missing = "??", callThresh = 0.9),
                n = 8)
  expect_is(compute_cf(marker_numeric, output = "markerCallfreq",
                       missing = NA_real_),
            class = "numeric")
})
