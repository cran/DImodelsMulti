library(DImodelsMulti)

test_that("Info Criteria Matches", {

  model <- DImulti(prop = c("G1", "G2", "L1", "L2"), y = "Y", eco_func = c("Var", "un"),
                   unit_IDs = "Plot", theta = c(0.5, 1, 1.2), DImodel = "AV",
                   method = "REML", data = dataBEL)


  expect_equal(class(AIC(model)),
                  "numeric")

  expect_equal(class(BIC(model)),
                  "numeric")

  expect_equal(class(AICc(model)),
                  "numeric")

  expect_equal(class(BICc(model)),
                  "numeric")


  expect_equal(signif(AIC(model), digits = 5),
              561.58)

  expect_equal(signif(BIC(model), digits = 5),
               610.25)

  expect_equal(signif(AICc(model), digits = 5),
               577.20)

  expect_equal(signif(BICc(model), digits = 5),
               640.82)

})
