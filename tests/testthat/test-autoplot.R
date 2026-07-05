test_that("autoplot rejects confMeta objects with differing k_studies", {
  est <- c(0.2, 0.5); se <- c(0.1, 0.2)
  cm1 <- confMeta(est, se, fun = p_edgington, input_p = "greater")
  cm2 <- confMeta(est, se, fun = p_edgington, input_p = "greater",
                  k_studies = 3)
  expect_error(autoplot(cm1, cm2, type = "forest"),
               "same")
})