context("CallFreq input")
test_that("'missing_value' has the same type as the elements of 'x'", {
  expect_error(compute_cf(marker_numeric, output = "marker_names",
                          missing_value = "??"))

})

test_that("'callThresh' is required when 'markerNames' are expected", {
  expect_error(compute_cf(marker_character, output = "marker_names",
                          missing_value = "??"))
})

context("CallFreq output")
test_that("output is correct", {
  expect_is(compute_cf(marker_character, output = "marker_names",
                       missing_value = "??", call_threshold = 0.9),
            class = "character")
  expect_is(compute_cf(marker_numeric, output = "marker_names",
                       missing_value = NA_real_, call_threshold = 0.9),
            class = "character")
  expect_length(compute_cf(marker_character, output = "marker_names",
                       missing_value = "??", call_threshold = 0.9),
                n = 8)
  expect_is(compute_cf(marker_numeric, output = "marker_callfreq",
                       missing_value = NA_real_),
            class = "numeric")
})
