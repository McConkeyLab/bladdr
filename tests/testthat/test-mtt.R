dummy_mtt <- function(){
  data.frame(
    condition = factor("a"),
    nm562 = c(0.177,  0.126,  0.0516, 0.0483, 0.0469, 0.0505),
    nm660 = c(0.0467, 0.0442, 0.0389, 0.0395, 0.0399, 0.041),
    drug  = c(1e-4,   1,      10,     100,    1e3,    1e4)
  )
}

dummy_fit <- function() {
  mtt_calc(dummy_mtt())
}

dummy_curve <- function() {
  mtt <- dummy_mtt()
  fit <- dummy_fit()
  make_curve(fit$fit[[1]], c(mtt$drug), length_out = 10)
}

# make_curve -------------------------------------------------------------------
test_that("Curve is length_out", {
  curve <- dummy_curve()
  expect_equal(nrow(curve), 10)
})

test_that("Curve starts and ends at correct x bounds", {
  fit <- dummy_fit()
  curve <- dummy_curve()
  expect_equal(min(fit$drug), curve[1,1])
  expect_equal(max(fit$drug), curve[nrow(curve), 1])
})

test_that("NULL fit produces NA y but not NA x", {
  mtt <- dummy_mtt()
  expect_true(
    all(make_curve(NULL, c(mtt$drug), length_out = 10) |>
          dplyr::pull(y) |> is.na())
  )
  expect_true(
    !any(make_curve(NULL, c(mtt$drug), length_out = 10) |>
           dplyr::pull(x) |> is.na())
  )
})

# get_ic -----------------------------------------------------------------------
test_that("IC50 is consistent", {
  fit <- dummy_fit()
  ic50 <- fit$fit[[1]] |> get_ic(ic_pct = 50)
  expect_equal(dplyr::pull(ic50, ic_value) |> round(), 1)
})

test_that("IC50 is greater than IC25 for typical curves", {
  fit <- dummy_fit()
  ic50 <- fit$fit[[1]] |> get_ic(ic_pct = 50)
  ic25 <- fit$fit[[1]] |> get_ic(ic_pct = 25)
  expect_gt(dplyr::pull(ic50, ic_value), dplyr::pull(ic25, ic_value))
})

test_that("NULL fit produces NA output", {
  expect_true(is.na(get_ic(NULL, ic_pct = 50)))
})


# Sanitize_drug_conc -----------------------------------------------------------
test_that("0 is converted to something small", {
  expect_equal(
    sanitize_drug_conc(c(0, 1, 10)),
    c(0.0001, 1, 10)
  )
  expect_equal(
    sanitize_drug_conc(c(0, 1, 2)),
    c(0.0625, 1, 2)
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
    c(1, 0.0625, 2)
  )
})

test_that("New 0 is correct even when numerator is not 1", {
  expect_equal(
    sanitize_drug_conc(c(0.01, 0.1, 0, 1, 10, 100)),
    c(0.01, 0.1, 1e-4, 1, 10, 100)
  )
})
