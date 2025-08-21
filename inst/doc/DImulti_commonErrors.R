## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

options(crayon.enabled = TRUE)

ansi_aware_handler <- function(x, options)
{
  paste0(
    "<pre class=\"r-output\"><code>",
    fansi::sgr_to_html(x = x, warn = FALSE, term.cap = "256"),
    "</code></pre>"
  )
}

old_hooks <- fansi::set_knit_hooks(knitr::knit_hooks,
                                 which = c("output", "message", "error", "warning"))
knitr::knit_hooks$set(
  output = ansi_aware_handler,
  message = ansi_aware_handler,
  warning = ansi_aware_handler,
  error = ansi_aware_handler
)

## ----setup--------------------------------------------------------------------
library(DImodelsMulti)
data("dataBEL"); data("dataSWE"); data("simMV"); data("simRM"); data("simMVRM")

## ----y_univar, eval = FALSE---------------------------------------------------
# DImulti(y = "YIELD", ..., data = dataSWE)
# DImulti(y = 9, ..., data = dataSWE)

## ----y_dataSWE, eval = TRUE, echo = FALSE-------------------------------------
dataSWE[1:6, "YIELD", drop = FALSE]

## ----y_MVstack, eval = FALSE--------------------------------------------------
# DImulti(y = "Y", ..., data = dataBEL)
# DImulti(y = 9, ..., data = dataBEL)

## ----y_dataBEL, eval = TRUE, echo = FALSE-------------------------------------
dataBEL[1:6, "Y", drop = FALSE]

## ----y_MVwide, eval = FALSE---------------------------------------------------
# DImulti(y = c("Y1", "Y2", "Y3", "Y4"), ..., data = simMV)
# DImulti(y = 9:12, ..., data = simMV)

## ----y_simMV, eval = TRUE, echo = FALSE---------------------------------------
simMV[1:6, c("Y1", "Y2", "Y3", "Y4"), drop = FALSE]

## ----ecoFunc_MVstack, eval = FALSE--------------------------------------------
# DImulti(eco_func = c("Var", "UN"), ..., data = dataBEL)

## ----ecoFunc_dataBEL, eval = TRUE, echo = FALSE-------------------------------
dataBEL[1:6, "Var", drop = FALSE]

## ----ecoFunc_MVwide, eval = FALSE---------------------------------------------
# DImulti(eco_func = c("NA", "UN"), ..., data = simMV)

## ----time, eval = FALSE-------------------------------------------------------
# DImulti(time = c("YEARN", "AR1"), ..., data = dataSWE)

## ----time_dataSWE, eval = TRUE, echo = FALSE----------------------------------
dataSWE[c(1, 49, 97, 2, 50, 98), "YEARN", drop = FALSE]

## ----unit_IDs, eval = FALSE---------------------------------------------------
# DImulti(unit_IDs = "plot", ..., data = simMVRM)
# DImulti(unit_IDs = 1, ..., data = simMVRM)

## ----unit_IDs_simVMRM, eval = TRUE, echo = FALSE------------------------------
simMVRM[1:6, "plot", drop = FALSE]

## ----prop, eval = FALSE-------------------------------------------------------
# DImulti(unit_IDs = paste0("p", 1:4), ..., data = simMVRM)
# DImulti(unit_IDs = 2:5, ..., data = simMVRM)

## ----prop_simVMRM, eval = TRUE, echo = FALSE----------------------------------
simMVRM[1:6, 2:5, drop = FALSE]

## ----data, eval = FALSE-------------------------------------------------------
# DImulti(..., data = simMVRM)
# DImulti(..., data = simMVRM)

## ----data_simVMRM, eval = TRUE, echo = FALSE----------------------------------
simMVRM[1:6, , drop = FALSE]

## ----DImodel, eval = FALSE----------------------------------------------------
# DImulti(DImodel = "FULL", ...)
# DImulti(DImodel = "AV", ...)

## ----FG, eval = FALSE---------------------------------------------------------
# DImulti(prop = c("G1", "G2", "L1", "L2"), DImodel = "FG",
#         FG = c("Grass", "Grass", "Legume", "Legume"), ..., data = dataBEL)

## ----ID, eval = FALSE---------------------------------------------------------
# DImulti(prop = c("G1", "G2", "L1", "L2"), ID = c("Grass", "Grass", "Legume", "Legume"),
#         ..., data = dataBEL)

## ----extra_fixed, eval = FALSE------------------------------------------------
# DImulti(extra_fixed = ~DENS, ..., data = dataSWE)
# DImulti(extra_fixed = ~DENS + TREAT, ..., data = dataSWE)
# DImulti(extra_fixed = ~1*DENS, ..., data = dataSWE)
# DImulti(extra_fixed = ~1:(DENS+TREAT), ..., data = dataSWE)
# DImulti(extra_fixed = ~1:DENS:TREAT, ..., data = dataSWE)

## ----extra_fixed_dataSWE, eval = TRUE, echo = FALSE---------------------------
dataSWE[1:6, c("DENS", "TREAT"), drop = FALSE]

## ----estimate_theta, eval = FALSE---------------------------------------------
# DImulti(estimate_theta = TRUE, ...)
# DImulti(estimate_theta = FALSE, ...)

## ----theta, eval = FALSE------------------------------------------------------
# DImulti(theta = 1, ...)
# DImulti(theta = c(0.5, 1, 1.2), ...)

## ----method, eval = FALSE-----------------------------------------------------
# DImulti(method = "REML", ...)
# DImulti(method = "ML", ...)

## ----singular_dataSWE, eval = TRUE, error = TRUE------------------------------
try({
DImulti(prop = 5:8, y = "YIELD", time = c("YEARN", "CS"), unit_IDs = "PLOT", data = dataSWE, 
        DImodel = "ID", extra_fixed = ~DENS+TREAT, method = "REML")
})

## ----singular_dataSWE_DENSTREAT, eval = TRUE, echo = FALSE--------------------
dataSWE[c(1, 16, 31, 40), c("DENS", "TREAT"), drop = FALSE]

## ----singular_dataSWE_fix, eval = TRUE----------------------------------------
DImulti(prop = 5:8, y = "YIELD", time = c("YEARN", "CS"), unit_IDs = "PLOT", data = dataSWE, 
        DImodel = "ID", extra_fixed = ~1:(DENS+TREAT), method = "REML")

## ----singular_dataBEL_Works, eval = TRUE--------------------------------------
model1 <- DImulti(prop = 2:5, y = "Y", eco_func = c("Var", "UN"), unit_IDs = 1, data = dataBEL,
                  FG = c("Grass", "Grass", "Leg", "Leg"), DImodel = "FG", extra_fixed = ~Density, 
                  method = "ML", theta = 1)

## ----comparison_dataBEL, eval = TRUE------------------------------------------
model1 <- DImulti(prop = 2:5, y = "Y", eco_func = c("Var", "UN"), unit_IDs = 1, data = dataBEL,
                  FG = c("Grass", "Grass", "Leg", "Leg"), DImodel = "FG", 
                  method = "REML")

model2 <- DImulti(prop = 2:5, y = "Y", eco_func = c("Var", "UN"), unit_IDs = 1, data = dataBEL,
                  FG = c("Grass", "Grass", "Leg", "Leg"), DImodel = "FG", extra_fixed = ~Density,
                  method = "REML")

anova(model1, model2)

