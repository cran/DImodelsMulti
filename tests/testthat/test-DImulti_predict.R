library(DImodelsMulti)

test_that("DImulti Prediction Tests", {

  model <- DImulti(y = 6:8, eco_func = c("Na", "un"), time = c("time", "CS"), unit_IDs = 1,
                   prop = 2:5, data = simMVRM, DImodel = "AV",
                   estimate_theta = TRUE, method = "REML")

  expect_s3_class(predict(model, stacked = FALSE),
                  "data.frame")

  expect_s3_class(predict(model, newdata = simMVRM[1, ], stacked = FALSE),
                  "data.frame")

  expect_s3_class(predict(model, stacked = TRUE),
                  "data.frame")

  expect_s3_class(predict(model, newdata = simMVRM[1:3, ], stacked = TRUE),
                  "data.frame")

  expect_s3_class(predict(model, newdata = simMVRM[1, ], stacked = TRUE),
                  "data.frame")


  expect_equal(colnames(predict(model, newdata = simMVRM[1:6, ])),
               c("plot", "Yvalue", "Ytype"))

  expect_equal(colnames(predict(model, newdata = simMVRM[1:6, ], stacked = FALSE)),
               c("plot", "Y1:1", "Y2:1", "Y3:1"))


  expect_warning(predict(model, newdata = simMVRM[1:6, -1]),
                 "The column containing unit_IDs has not been supplied through newdata. This column is required as a grouping factor for the covarying responses, although its value does not matter as there is no between subject effect included. Defaulting to row numbers.")


  expect_error(predict(model, newdata = simMVRM[1:6, -2:-3]),
               "The following initial proportion columns are missing from 'newdata':\np1, p2")

})
