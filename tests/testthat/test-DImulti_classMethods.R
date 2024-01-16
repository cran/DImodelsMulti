library(DImodelsMulti)

test_that("S3 Methods Checks", {

  model <- DImulti(y = 6:8, eco_func = c("Na", "un"), time = c("time", "CS"), unit_IDs = 1,
                   prop = 2:5, data = simMVRM, DImodel = "AV",
                   estimate_theta = TRUE, method = "REML")

  expect_s3_class(print(model),
                  "DImulti")

  expect_s3_class(summary(model),
                  "summary.DImulti")

  expect_s3_class(print(summary(model)),
                  "summary.DImulti")


  expect_match(toString(capture.output(print(model))),
               "^Note.*Generalized.*Table.*Degrees.*")

  expect_match(toString(capture.output(print(summary(model)))),
               "^Generalized.*Model.*Correlation.*Coefficients.*Theta.*Correlation.*Standardized.*")

})
