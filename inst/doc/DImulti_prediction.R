## ---- include = FALSE---------------------------------------------------------
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

## ----DImulti_modelEx----------------------------------------------------------
modelFinal <- DImulti(y = c("Y1", "Y2", "Y3"), eco_func = c("NA", "UN"), time = c("time", "CS"),
                    unit_IDs = 1, prop = 2:5, data = simMVRM, DImodel = "AV", method = "REML",
                    estimate_theta = TRUE)
print(modelFinal)

## ----predict_layout, eval=FALSE-----------------------------------------------
#  predict.DImulti(object, newdata = NULL, stacked = TRUE, ...)

## ----predict_default----------------------------------------------------------
head(predict(modelFinal))

## ----predict_wide-------------------------------------------------------------
head(predict(modelFinal, stacked = FALSE))

## ----predict_subset-----------------------------------------------------------
predict(modelFinal, newdata = simMVRM[c(1, 4, 7, 10, 21), ])

## ----predict_newSim-----------------------------------------------------------
newSim <- data.frame(plot = c(1, 2),
                     p1 = c(0.25, 0.6),
                     p2 = c(0.25, 0.2),
                     p3 = c(0.25, 0.1),
                     p4 = c(0.25, 0.1)) 

predict(modelFinal, newdata = newSim)

## ----predict_Y1---------------------------------------------------------------
newSim <- data.frame(plot = c(1, 2),
                     p1 = c(0.25, 0.6),
                     p2 = c(0.25, 0.2),
                     p3 = c(0.25, 0.1),
                     p4 = c(0.25, 0.1),
                     Y1 = 0) 

predict(modelFinal, newdata = newSim)

## ----predict_newSim_missingID-------------------------------------------------
newSim <- data.frame(p1 = c(0.25, 0.6),
                     p2 = c(0.25, 0.2),
                     p3 = c(0.25, 0.1),
                     p4 = c(0.25, 0.1)) 

predict(modelFinal, newdata = newSim)

## ----predict_newSim_merge-----------------------------------------------------
newSim <- data.frame(plot = c(1, 2),
                     p1 = c(0.25, 0.6),
                     p2 = c(0.25, 0.2),
                     p3 = c(0.25, 0.1),
                     p4 = c(0.25, 0.1)) 

preds <- predict(modelFinal, newdata = newSim, stacked = FALSE)

merge(newSim, preds, by = "plot")

## ----predict_newSim_aggregate-------------------------------------------------
newSim <- data.frame(plot = c(1, 1),
                     p1 = c(0.25, 0.6),
                     p2 = c(0.25, 0.2),
                     p3 = c(0.25, 0.1),
                     p4 = c(0.25, 0.1)) 

predict(modelFinal, newdata = newSim, stacked = FALSE)

