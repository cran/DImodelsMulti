library(DImodelsMulti)

test_that("Missing/incorrect parameters", {
  expect_error(DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
                       unit_IDs = 2, DImodel = "ID", ID = c("G", "G", "L", "L"),
                       method = "REML"),
               "You must supply a dataframe through the argument 'data'.\n")

  expect_error(DImulti(y = c("YIELD"), time = c("YEARN", "CS"),
                       unit_IDs = 2, DImodel = "ID", ID = c("G", "G", "L", "L"),
                       method = "REML", data = dataSWE),
               "You must supply species proportion variable names or column indices through the argument 'prop'.\n")

  expect_error(DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
                       unit_IDs = 2, DImodel = "ID", ID = c("G", "G", "L", "L"),
                       method = "REML", data = cbind(dataSWE, E = 0)),
               "You must not have any columns named \"AV\" or \"E\", \"NA\" or \"NULL\" in your dataset provided through the parameter 'data'")

  expect_error(DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
                       unit_IDs = 2, DImodel = "ID", ID = c("G", "G", "L", "L"), FG = c("Grass", "Legume"),
                       method = "REML", data = dataSWE),
               "The number of functional groups supplied does not match the number of species supplied through 'prop'")

  expect_error(DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
                       unit_IDs = 2, DImodel = "ID", ID = c("G", "G", "L", "L"),
                       extra_fixed = "a string", method = "REML", data = dataSWE),
               "You must supply a formula through the argument 'extra_fixed'\n")

  expect_error(DImulti(prop = 5:8, time = c("YEARN", "CS"),
                       unit_IDs = 2, DImodel = "ID", ID = c("G", "G", "L", "L"),
                       method = "REML", data = dataSWE),
               "You must supply response variable names or column indices through the argument 'y'.\n")

  expect_error(DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
                       DImodel = "ID", ID = c("G", "G", "L", "L"),
                       method = "REML", data = dataSWE),
               "You must either supply an ID variable name or column index through the argument 'unit_IDs'.\n")

  expect_error(DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
                       unit_IDs = 2, DImodel = "ID", ID = c("G", "G", "L", "L"),
                       method = "String", data = dataSWE),
               "You must either enter \"REML\" or \"ML\" as a string for the parameter 'method'")


  expect_s3_class(DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
                    unit_IDs = 2, DImodel = "ID", ID = c("G", "G", "L", "L"),
                    method = "REML", data = dataSWE),
            "DImulti")

  expect_s3_class(DImulti(prop = c("G1", "G2", "L1", "L2"), y = "Y", eco_func = c("Var", "un"),
                    unit_IDs = "Plot", theta = c(0.5, 1, 1.2), DImodel = "FG",
                    FG = c("Grass", "Grass", "Legume", "Legume"), extra_fixed = ~ Density,
                    method = "REML", data = dataBEL),
            "DImulti")

  expect_s3_class(DImulti(y = 6:8, eco_func = c("Na", "un"), time = c("time", "CS"), unit_IDs = 1,
                    prop = 2:5, data = simMVRM, DImodel = "AV",
                    estimate_theta = TRUE, method = "REML"),
            "DImulti")

  expect_s3_class(DImulti(y = 6:8, eco_func = c("Na", "un"), time = c("time", "CS"), unit_IDs = 1,
                          prop = 2:5, data = simMVRM, DImodel = "AV",
                          estimate_theta = TRUE, method = "REML")$corr$`Multivariate`,
                  "corSymm")

  expect_s3_class(DImulti(y = 6:8, eco_func = c("Na", "un"), time = c("time", "CS"), unit_IDs = 1,
                          prop = 2:5, data = simMVRM, DImodel = "AV",
                          estimate_theta = TRUE, method = "REML")$corr$`Repeated Measure`,
                  "corCompSymm")


  expect_equal(signif(DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
                           unit_IDs = 2, DImodel = "AV", ID = c("G", "G", "L", "L"),
                           estimate_theta = TRUE, method = "REML", data = dataSWE)$theta, digits = 2),
               6.8e-05)

  expect_equal(DImulti(prop = 5:8, y = c("YIELD"), time = c("YEARN", "CS"),
                              unit_IDs = 2, DImodel = "AV", ID = c("G", "G", "L", "L"),
                              theta = 1.5, method = "REML", data = dataSWE)$theta,
               1.5)

  expect_equal(DImulti(prop = c("G1", "G2", "L1", "L2"), y = "Y", eco_func = c("Var", "un"),
                       unit_IDs = "Plot", theta = c(0.5, 1, 1.2), DImodel = "AV",
                       method = "REML", data = dataBEL)$theta,
               c("N" = 0.5, "Sown" = 1, "Weed" = 1.2))
})


