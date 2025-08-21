## ----include=FALSE------------------------------------------------------------
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

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(DImodelsMulti)
library(dplyr)
library(nlme)
library(ggplot2)

## ----data_head----------------------------------------------------------------
head(simMVRM)

## ----data_group---------------------------------------------------------------
simMVRM_group <-  dplyr::summarise(dplyr::group_by(simMVRM, time),
                                   Y1 = mean(Y1),
                                   Y2 = mean(Y2),
                                   Y3 = mean(Y3),
                                   MFindex = mean(Y1 + Y2 + Y3))
simMVRM_group

## ----data_hist, out.width="75%", fig.alt="A histogram of the raw data"--------
hist(simMVRM[which(simMVRM$time == 1), ]$Y1, main = "Y1 at time 1", xlab = "Y1")

## ----DImulti_modelSTR---------------------------------------------------------
modelSTR <- DImulti(y = c("Y1", "Y2", "Y3"), eco_func = c("NA", "UN"), time = c("time", "CS"),
                    unit_IDs = 1, prop = 2:5, data = simMVRM, DImodel = "STR", method = "ML")

## ----DImulti_modelSTR_print---------------------------------------------------
print(modelSTR)

## ----DImulti_modelID----------------------------------------------------------
modelID <- DImulti(y = c("Y1", "Y2", "Y3"), eco_func = c("NA", "UN"), time = c("time", "CS"),
                    unit_IDs = 1, prop = 2:5, data = simMVRM, DImodel = "ID", method = "ML")

## ----DImulti_modelID_tTable---------------------------------------------------
summary(modelID)$tTable

## ----DImulti_STR_ID_LRT-------------------------------------------------------
anova(modelSTR, modelID)

## ----DImulti_modelAV----------------------------------------------------------
modelAV <- DImulti(y = c("Y1", "Y2", "Y3"), eco_func = c("NA", "UN"), time = c("time", "CS"),
                    unit_IDs = 1, prop = 2:5, data = simMVRM, DImodel = "AV", method = "ML")

coef(modelAV)
anova(modelID, modelAV)

## ----DImulti_modelAV_theta----------------------------------------------------
modelAV_theta <- DImulti(y = c("Y1", "Y2", "Y3"), eco_func = c("NA", "UN"), time = c("time", "CS"),
                    unit_IDs = 1, prop = 2:5, data = simMVRM, DImodel = "AV", method = "ML",
                    estimate_theta = TRUE)

thetaVals <- modelAV_theta$theta
thetaVals

AICc(modelAV) 
AICc(modelAV_theta)

## ----DImulti_modelADD---------------------------------------------------------
modelADD <- DImulti(y = c("Y1", "Y2", "Y3"), eco_func = c("NA", "UN"), time = c("time", "CS"),
                    unit_IDs = 1, prop = 2:5, data = simMVRM, DImodel = "ADD", method = "ML",
                    theta = thetaVals)

modelADD$coefficients
anova(modelAV_theta, modelADD)

## ----DImulti_modelADD_treat, eval=FALSE---------------------------------------
# modelAV_treat1 <- DImulti(y = c("Y1", "Y2", "Y3"), eco_func = c("NA", "UN"), time = c("time", "CS"),
#                     unit_IDs = 1, prop = 2:5, data = simMVRM, DImodel = "ADD", method = "ML",
#                     theta = thetaVals, extra_fixed = ~treat)
# 
# modelAV_treat2 <- DImulti(y = c("Y1", "Y2", "Y3"), eco_func = c("NA", "UN"), time = c("time", "CS"),
#                     unit_IDs = 1, prop = 2:5, data = simMVRM, DImodel = "ADD", method = "ML",
#                     theta = thetaVals, extra_fixed = ~1:treat)

## ----DImulti_modelAV_ID-------------------------------------------------------
modelAV_ID <- DImulti(y = c("Y1", "Y2", "Y3"), eco_func = c("NA", "UN"), time = c("time", "CS"),
                    unit_IDs = 1, prop = 2:5, data = simMVRM, DImodel = "AV", method = "ML",
                    theta = thetaVals, ID = c("Group1", "Group1", "Group2", "Group2"))

summary(modelAV_ID)$tTable
anova(modelAV_ID, modelAV)

## ----DImulti_modelFinal-------------------------------------------------------
modelFinal <- DImulti(y = c("Y1", "Y2", "Y3"), eco_func = c("NA", "UN"), time = c("time", "CS"),
                    unit_IDs = 1, prop = 2:5, data = simMVRM, DImodel = "AV", method = "REML",
                    theta = thetaVals)

summary(modelFinal)

## ----DImulti_modelFinal_tTable------------------------------------------------
summary(modelFinal)$tTable

## ----DImulti_modelFinal_varCovs, eval=FALSE-----------------------------------
# nlme::getVarCov(modelFinal)
# 
# modelFinal$vcov

## ----DImulti_modelFinal_getVarCov, echo=FALSE---------------------------------
nlme::getVarCov(modelFinal)

## ----DImulti_modelFinal_predict-----------------------------------------------
comms <- simMVRM[c(1, 4, 7, 10, 21), ]
print(comms)


commPred <- predict(modelFinal, newdata = comms)
commPred

## ----DImulti_modelFinal_grouped, out.width="75%", fig.alt="An example histogram"----
ggplot2::ggplot(commPred, ggplot2::aes(fill=Ytype, y=Yvalue, x=plot)) +
  ggplot2::geom_bar(position="dodge", stat="identity")

