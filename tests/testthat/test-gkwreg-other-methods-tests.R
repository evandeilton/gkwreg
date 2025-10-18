library(testthat)
library(gkwreg)

# Test data setup
data("GasolineYield")

test_that("Basic model update works for simple formulas", {
  # Fit initial model
  m0 <- gkwreg(yield ~ 1, data = GasolineYield, family = "kw")

  # Add predictor
  m1 <- update(m0, . ~ . + temp)
  expect_s3_class(m1, "gkwreg")
  expect_gt(AIC(m0), AIC(m1)) # Better fit expected

  # Remove predictor (back to intercept-only)
  m2 <- update(m1, . ~ . - temp)
  expect_equal(AIC(m0), AIC(m2), tolerance = 1e-6)
})

test_that("Multi-part formula updates work correctly", {
  # Start with two-part model
  m0 <- gkwreg(yield ~ 1 | 1, data = GasolineYield, family = "kw")

  # Add temp to alpha only
  m1 <- update(m0, . ~ . + temp | .)
  expect_s3_class(m1, "gkwreg")

  # Add batch to beta only
  m2 <- update(m1, . ~ . | . + batch)
  expect_s3_class(m2, "gkwreg")

  # Verify both predictors are included appropriately
  expect_true(any(grepl("temp", names(coef(m1)))))
  expect_true(any(grepl("batch", names(coef(m2)))))
})

test_that("Three-part formula updates function properly", {
  # Initial three-part model
  m0 <- gkwreg(yield ~ 1 | 1 | 1, data = GasolineYield, family = "gkw")

  # Update each part separately
  m1 <- update(m0, . ~ . + temp | . | .) # alpha
  m2 <- update(m1, . ~ . | . + batch | .) # beta
  m3 <- update(m2, . ~ . | . | . + temp) # gamma

  # Check that all parts were updated correctly
  expect_s3_class(m3, "gkwreg")
  formula_str <- as.character(formula(m3))
  expect_true(grepl("temp", formula_str[3])) # alpha part
  expect_true(grepl("batch", formula_str[3])) # beta part
  expect_true(grepl("temp", formula_str[3])) # gamma part
})

test_that("Formula extraction works for different formula types", {
  # Simple formula
  fit1 <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")
  expect_s3_class(formula(fit1), "formula")

  # Two-part formula
  fit2 <- gkwreg(yield ~ temp | batch, data = GasolineYield, family = "kw")
  expect_s3_class(formula(fit2), "Formula")

  # Five-part formula
  fit3 <- gkwreg(yield ~ temp | batch | temp | 1 | 1,
    data = GasolineYield, family = "gkw"
  )
  expect_s3_class(formula(fit3), "Formula")
})

test_that("Model components can be extracted properly", {
  fit <- gkwreg(yield ~ temp + batch, data = GasolineYield, family = "kw")

  # Test model frame extraction
  mf <- model.frame(fit)
  expect_s3_class(mf, "data.frame")
  expect_true(all(c("yield", "temp", "batch") %in% names(mf)))

  # Test model matrix extraction
  mm <- model.matrix(fit)
  expect_type(mm, "double")
  expect_true(is.matrix(mm))

  # Test terms extraction
  tt <- terms(fit)
  expect_s3_class(tt, "terms")

  # Test response extraction
  y <- response(fit)
  expect_type(y, "double")
  expect_length(y, nrow(GasolineYield))

  # Test family extraction
  fam <- family(fit)
  expect_type(fam, "character")
  expect_true(fam %in% c("kw", "gkw"))
})

test_that("Update with new data works correctly", {
  # Create subset of data
  data_subset <- GasolineYield[1:20, ]

  # Fit on full data
  fit1 <- gkwreg(yield ~ temp, data = GasolineYield, family = "kw")

  # Update with new data
  fit2 <- update(fit1, data. = data_subset)

  # Should have different number of observations
  expect_equal(nobs(fit2), 20)
  expect_equal(nobs(fit1), nrow(GasolineYield))
})
