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
library(DImodels)

## ----thetaOptions, eval = FALSE-----------------------------------------------
#  DImulti(..., theta = c(0.5, 1, 1.2))
#  DImulti(..., estimate_theta = TRUE)
#  DImulti(..., theta = 1, estimate_theta = FALSE)

## ----DI_data_value, eval = FALSE----------------------------------------------
#  DI_data(..., theta = 1.2, what = "FULL")

## ----DI_estimate, eval = FALSE------------------------------------------------
#  DI(..., estimate_theta = TRUE, DImodel = "FULL")$coefficients[["theta"]]

## ----optimTheta, eval = FALSE-------------------------------------------------
#  f <- function(theta_val)
#  {
#    fit <- DImulti(prop = 2:5, y = 6:8, eco_func = c("na", "un"), time = c("time", "ar1"),
#                   unit_IDs = 1, theta = theta_val, DImodel = "FULL", method = "ML",
#                   data = simMVRM)
#  
#  	return(-as.numeric(logLik(fit)))
#  }
#  	
#  optim(c(1,1,1), f, hessian = FALSE, lower = c(.1,.1,.1), upper = c(1.5,1.5,1.5),
#        method = "L-BFGS-B")

