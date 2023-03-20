test_that("0 is converted to something small", {
  expect_equal(
    sanitize_drug_conc(c(0, 1, 10)),
    c(0.001, 1, 10)
  )
  expect_equal(
    sanitize_drug_conc(c(0, 1, 2)),
    c(0.125, 1, 2)
    )
})

test_that("No 0 means no conversion", {
  expect_equal(
    sanitize_drug_conc(c(0.1, 1, 2)),
    c(0.1, 1, 2)
  )
})

test_that("Drug concentration is not reordered", {
  expect_equal(
    sanitize_drug_conc(c(1, 0, 2)),
    c(1, 0.125, 2)
  )
})
