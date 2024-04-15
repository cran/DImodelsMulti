###########################################################################################################################################################################################
#Internal function to set correlation structure (corr variable) for use in DImulti_fit
#' @keywords internal
#' @noRd

DI_vcov <- function(model, MVflag, Timeflag, funcCorr = "NA", timeCorr = "NA", formulaStart, formulaEnd, prop, ID, Yfunc, Yvalue, unitIDs, method, timeCol = "NA", data, theta, exFixflag)
{
  corrMV <- corrT <- MVvcov <- Tvcov <- NULL

  #For MV or Time, not both (simpler -> faster)
  if((MVflag & !Timeflag) | (!MVflag & Timeflag))
  {
    corrform <- if(MVflag) funcCorr else timeCorr

    redo <- FALSE
    repeat
    {
      corr <- switch(
        tolower(corrform),
        "un" = {
          redo <- FALSE
          if(MVflag)
          {
            corr <- corrMV <- nlme::corSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))
          }
          else
          {
            corr <- corrT <- nlme::corSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))
          }
        },
        "cs" = {
          redo <- FALSE
          if(MVflag)
          {
            corr <- corrMV <- nlme::corCompSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))
          }
          else
          {
            corr <- corrT <- nlme::corCompSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))
          }
        },
        "ar1" = {
          if(MVflag)
          {
            warning(crayon::bold(crayon::red("AR(1) is a nonsensical structure for ecosystem functions, defaulting to \"un\"")))
            redo <- TRUE
            corrform <- "un"
            corr <- NA
          }
          else
          {
            corr <- corrT <- nlme::corAR1(form = stats::as.formula(paste("~ 0 + as.numeric(", timeCol, ") | ", paste(unitIDs, collapse = "/"))))
          }
        },
        "na" = {
          warning(crayon::bold(crayon::red("You have not supplied a correlation structure through the argument \"ecoFunc\", defaulting to
                                           \"un\"")))
          redo <- TRUE
          corrform <- "un"
          corr <- NA
        },
        {
          warning(crayon::bold(crayon::red("Chosen correlation structure not available, defaulting to \"un\"")))
          redo <- TRUE
          corrform <- "un"
          corr <- NA
        }
      )

      if(!redo) { break }
    }

    #END of if
    return (list(corr, corrform, corrMV, corrT))
    ##########


  }
  else #Both MV & Rep
  {
    corrform <- paste0(funcCorr, "@", timeCorr)

    redo <- FALSE
    repeat
    {
      corr <- switch(
        tolower(corrform),
        "un@un" = {
          #Usual UN structure for both
          un <- nlme::corSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))

          #Fit model here for func and time, holding other constant (saving second value as first is formula)
          models <- DImulti_vcovModels(model, MVflag, funcCorr = list(un, "un"), timeCorr = list(un, "un"), formulaStart, formulaEnd, unitIDs,
                                       prop, ID, Yfunc, Yvalue, method, timeCol, data, theta, exFixflag)
          MVmodel <- models[[1]][[1]]
          Tmodel <- models[[2]][[1]]

          # MVval <- coef(MVmodel$modelStruct$corStruct, unconstrained = FALSE)
          # Tval <- coef(Tmodel$modelStruct$corStruct, unconstrained = FALSE)
          #
          # #Initialize for func (hold time constant)
          # unMV <- nlme::corSymm(value = MVval, form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))), fixed = TRUE)
          # corrMV <- nlme::Initialize(unMV, data[which(data[, timeCol] == data[1, timeCol]), ])
          # unMV <- nlme::corMatrix(corrMV, corr = TRUE) #list of block diag
          #
          # #Initialize for time (hold func constant)
          # unT <- nlme::corSymm(value = Tval, form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))), fixed = TRUE)
          # corrT <- nlme::Initialize(unT, data[which(data[, Yfunc] == data[1, Yfunc]), ])
          # unT <- nlme::corMatrix(corrT, corr = TRUE) # list of block diag
          #
          # #sd for both MV and time
          # sdMV <- c()
          # for(i in unique(data[, Yfunc]))
          # {
          #   sdMV <- append(sdMV, sd(data[which(data[Yfunc] == i), "value"]))
          # }
          # sdT <- c()
          # for(i in unique(as.factor(data[, timeCol])))
          # {
          #   sdT <- append(sdT, sd(data[which(data[timeCol] == i), "value"]))
          # }
          #
          # #convert cors to covs
          # unMV <- MBESS::cor2cov(unMV[[1]], sd = sdMV)
          # unT <- MBESS::cor2cov(unT[[1]], sd = sdT)


          #Extract individual corStructs
          corrMV <- MVmodel$modelStruct$corStruct
          MVvcov <- nlme::getVarCov(MVmodel)

          corrT <- Tmodel$modelStruct$corStruct
          Tvcov <- nlme::getVarCov(Tmodel)


          unMV <- nlme::getVarCov(MVmodel)
          unT <- nlme::getVarCov(Tmodel)

          un <- kronecker(unMV, unT)

          #convert back to cor and fit into general structure
          un <- Matrix::cov2cor(un)

          #just in case, make positive-definite
          un <- Matrix::nearPD(un)$mat

          #Symmetrical on diag, so keep lower triangle
          un <- un[lower.tri(un)]

          #Use as fixed values for general structure (get corStruct object)
          un <- nlme::corSymm(value = un, form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))), fixed = TRUE)

          redo <- FALSE
          corr <- un
        },
        "un@cs" = {
          #Usual un and cs for both
          un <- nlme::corSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))
          cs <- nlme::corCompSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))

          #fit model here for func and time, holding other constant (saving second value as first is formula)
          models <- DImulti_vcovModels(model, MVflag, funcCorr = list(un, "un"), timeCorr = list(cs, "cs"), formulaStart, formulaEnd, unitIDs,
                                       prop, ID, Yfunc, Yvalue, method, timeCol, data, theta, exFixflag)
          MVmodel <- models[[1]][[1]]
          Tmodel <- models[[2]][[1]]

          corrMV <- MVmodel$modelStruct$corStruct
          MVvcov <- nlme::getVarCov(MVmodel)

          corrT <- Tmodel$modelStruct$corStruct
          Tvcov <- nlme::getVarCov(Tmodel)

          un <- nlme::getVarCov(MVmodel)
          cs <- nlme::getVarCov(Tmodel)

          uncs <- kronecker(un, cs)

          #convert back to cor and fit into general structure
          uncs <- Matrix::cov2cor(uncs)

          #just in case, make positive-definite
          uncs <- Matrix::nearPD(uncs)$mat

          #Symmetrical on diag, so keep lower triangle
          uncs <- uncs[lower.tri(uncs)]

          #Use as fixed values for general structure (get corStruct object)
          uncs <- nlme::corSymm(value = uncs, form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))), fixed = TRUE)

          redo <- FALSE
          corr <- uncs
        },

        "un@ar1" = {
          #Usual un and ar1 for both
          un <- nlme::corSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))
          ar1 <- nlme::corAR1(form = stats::as.formula(paste("~ 0 + as.numeric(", timeCol, ") | ", paste(unitIDs, collapse = "/"), "/", Yfunc)))

          #fit model here for func and time, holding other constant (saving second value as first is formula)
          models <- DImulti_vcovModels(model, MVflag, funcCorr = list(un, "un"), timeCorr = list(ar1, "ar1"), formulaStart, formulaEnd,
                                       unitIDs, prop, ID, Yfunc, Yvalue, method, timeCol, data, theta, exFixflag)
          MVmodel <- models[[1]][[1]]
          Tmodel <- models[[2]][[1]]

          corrMV <- MVmodel$modelStruct$corStruct
          MVvcov <- nlme::getVarCov(MVmodel)

          corrT <- Tmodel$modelStruct$corStruct
          Tvcov <- nlme::getVarCov(Tmodel)

          un <- nlme::getVarCov(MVmodel)
          ar1 <- nlme::getVarCov(Tmodel)

          unar1 <- kronecker(un, ar1)

          #convert back to cor and fit into general structure
          unar1 <- Matrix::cov2cor(unar1)

          #just in case, make positive-definite
          unar1 <- Matrix::nearPD(unar1)$mat

          #Symmetrical on diag, so keep lower triangle
          unar1 <- unar1[lower.tri(unar1)]

          #Use as fixed values for general structure (get corStruct object)
          unar1 <- nlme::corSymm(value = unar1, form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))), fixed = TRUE)

          redo <- FALSE
          corr <- unar1
        },
        "cs@cs" = {
          #usual cs for both
          cs <- nlme::corCompSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))

          #fit model here for func and time, holding other constant (saving second value as first is formula)
          models <- DImulti_vcovModels(model, MVflag, funcCorr = list(cs, "cs"), timeCorr = list(cs, "cs"), formulaStart, formulaEnd,
                                       unitIDs, prop, ID, Yfunc, Yvalue, method, timeCol, data, theta, exFixflag)
          MVmodel <- models[[1]][[1]]
          Tmodel <- models[[2]][[1]]

          corrMV <- MVmodel$modelStruct$corStruct
          MVvcov <- nlme::getVarCov(MVmodel)

          corrT <- Tmodel$modelStruct$corStruct
          Tvcov <- nlme::getVarCov(Tmodel)

          csMV <- nlme::getVarCov(MVmodel)
          csT <- nlme::getVarCov(Tmodel)

          cs <- kronecker(csMV, csT)

          #convert back to cor and fit into general structure
          cs <- Matrix::cov2cor(cs)

          #just in case, make positive-definite
          cs <- Matrix::nearPD(cs)$mat

          #Symmetrical on diag, so keep lower triangle
          cs <- cs[lower.tri(cs)]

          #Use as fixed values for general structure (get corStruct object)
          cs <- nlme::corSymm(value = cs, form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))), fixed = TRUE)

          redo <- FALSE
          corr <- cs
        },
        "cs@un" = {
          #Usual cs and un for both
          cs <- nlme::corCompSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))
          un <- nlme::corSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))

          #fit model here for func and time, holding other constant (saving second value as first is formula)
          models <- DImulti_vcovModels(model, MVflag, funcCorr = list(cs, "cs"), timeCorr = list(un, "un"), formulaStart, formulaEnd,
                                       unitIDs, prop, ID, Yfunc, Yvalue, method, timeCol, data, theta, exFixflag)
          MVmodel <- models[[1]][[1]]
          Tmodel <- models[[2]][[1]]

          corrMV <- MVmodel$modelStruct$corStruct
          MVvcov <- nlme::getVarCov(MVmodel)

          corrT <- Tmodel$modelStruct$corStruct
          Tvcov <- nlme::getVarCov(Tmodel)

          cs <- nlme::getVarCov(MVmodel)
          un <- nlme::getVarCov(Tmodel)

          csun <- kronecker(cs, un)

          #convert back to cor and fit into general structure
          csun <- Matrix::cov2cor(csun)

          #just in case, make positive-definite
          csun <- Matrix::nearPD(csun)$mat

          #Symmetrical on diag, so keep lower triangle
          csun <- csun[lower.tri(csun)]

          #Use as fixed values for general structure (get corStruct object)
          csun <- nlme::corSymm(value = csun, form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))), fixed = TRUE)

          redo <- FALSE
          corr <- csun
        },
        "cs@ar1" = {
          cs <- nlme::corCompSymm(form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))))
          ar1 <- nlme::corAR1(form = stats::as.formula(paste("~ 0 + as.numeric(", timeCol, ") | ", paste(unitIDs, collapse = "/"), "/", Yfunc)))

          #fit model here for func and time, holding other constant (saving second value as first is formula)
          models <- DImulti_vcovModels(model, MVflag, funcCorr = list(cs, "cs"), timeCorr = list(ar1, "ar1"), formulaStart, formulaEnd,
                                       unitIDs, prop, ID, Yfunc, Yvalue, method, timeCol, data, theta, exFixflag)
          MVmodel <- models[[1]][[1]]
          Tmodel <- models[[2]][[1]]

          corrMV <- MVmodel$modelStruct$corStruct
          MVvcov <- nlme::getVarCov(MVmodel)

          corrT <- Tmodel$modelStruct$corStruct
          Tvcov <- nlme::getVarCov(Tmodel)

          cs <- nlme::getVarCov(MVmodel)
          ar1 <- nlme::getVarCov(Tmodel)

          csar1 <- kronecker(cs, ar1)

          #convert back to cor and fit into general structure
          csar1 <- Matrix::cov2cor(csar1)

          #just in case, make positive-definite
          csar1 <- Matrix::nearPD(csar1)$mat

          #Symmetrical on diag, so keep lower triangle
          csar1 <- csar1[lower.tri(csar1)]

          #Use as fixed values for general structure (get corStruct object)
          csar1 <- nlme::corSymm(value = csar1, form = stats::as.formula(paste("~ 0 | ", paste(unitIDs, collapse = "/"))), fixed = TRUE)

          redo <- FALSE
          corr <- csar1

        },
        "ar1@un" = {
          warning(crayon::bold(crayon::red("AR(1) is a nonsensicle structure for ecoFunc, defaulting to \"un@un\"")))
          redo <- TRUE
          corrform <- "un@un"
        },
        "ar1@cs" = {
          warning(crayon::bold(crayon::red("AR(1) is a nonsensicle structure for ecoFunc, defaulting to \"un@cs\"")))
          redo <- TRUE
          corrform <- "un@cs"
        },
        "ar1@ar1" = {
          warning(crayon::bold(crayon::red("AR(1) is a nonsensicle structure for ecoFunc, defaulting to \"un@ar1\"")))
          redo <- TRUE
          corrform <- "un@ar1"
        },
        "na@un" = {
          warning(crayon::bold(crayon::red("You have not supplied a correlation structure through the argument \"ecoFunc\", defaulting to
                                           \"un\"")))
          redo <- TRUE
          corrform <- "un@un"
        },
        "na@cs" = {
          warning(crayon::bold(crayon::red("You have not supplied a correlation structure through the argument \"ecoFunc\", defaulting to
                                           \"un\"")))
          redo <- TRUE
          corrform <- "un@cs"
        },
        "na@ar1" = {
          warning(crayon::bold(crayon::red("You have not supplied a correlation structure through the argument \"ecoFunc\", defaulting to
                                           \"un\"")))
          redo <- TRUE
          corrform <- "un@ar1"
        },
        "na@na" = {
          warning(crayon::bold(crayon::red("You have not supplied a correlation structure through the argument \"ecoFunc\" or \"time\",
                                           defaulting to \"un@un\"")))
          redo <- TRUE
          corrform <- "un@un"
        },
        "un@na" = {
          warning(crayon::bold(crayon::red("You have not supplied a correlation structure through the argument \"time\", defaulting to
                                           \"un\"")))
          redo <- TRUE
          corrform <- "un@un"
        },
        "cs@na" = {
          warning(crayon::bold(crayon::red("You have not supplied a correlation structure through the argument \"time\", defaulting to
                                           \"un\"")))
          redo <- TRUE
          corrform <- "cs@un"
        },
        {
          warning(crayon::bold(crayon::red("Chosen correlation structure not available, defaulting to \"un@un\"")))
          redo <- TRUE
          corrform <- "un@un"
        }
      )

      if(!redo) { break }
    }

    #END of else
    return (list(corr, corrform, corrMV, corrT, MVvcov, Tvcov))
    ############
  }

  #END of method
  return (list(corr, corrform, corrMV, corrT, MVvcov, Tvcov))
  ##############
}


#To fit preemptively models using either Time or MV (when both given)
DImulti_vcovModels <- function(model, MVflag, funcCorr, timeCorr, formulaStart, formulaEnd, unitIDs, prop, ID, Yfunc, Yvalue, method,
                               timeCol, data, theta, exFixflag)
{
  if(tolower(funcCorr[[2]]) == "un")
  {
    weightMV <- nlme::varIdent(form = stats::as.formula(paste0("~ 0 |", Yfunc)))
  }
  else
  {
    weightMV <- NULL
  }
  if(tolower(timeCorr[[2]]) == "un")
  {
    weightT <- nlme::varIdent(form = stats::as.formula(paste0("~ 0 | ", timeCol)))
  }
  else
  {
    weightT <- NULL
  }

  funcCorr <- funcCorr[[1]]
  timeCorr <- timeCorr[[1]]

  dataMV <- data[which(data[, timeCol] == data[1, timeCol]), ]
  dataT <- data[which(data[, Yfunc] == data[1, Yfunc]), ]

  dataMV <- dataMV[order(dataMV[[unitIDs]], dataMV[[Yfunc]]), ]
  dataT <- dataT[order(dataT[[unitIDs]]), ]

  formulaStartT <- paste0(Yvalue, "~ 0 + ", timeCol, ":((")
  #remove any mention of func
  formulaEndT <- formulaEnd
#  formulaEndT <- gsub(paste0(Yfunc, "\\s*[:|\\*]"), "", formulaEnd)
#  formulaEndT <- gsub(paste0("[:|\\*]\\s*", Yfunc), "", formulaEndT)
#  formulaEndT <- gsub(paste0("\\+\\s*", Yfunc), "", formulaEndT)

  formulaStartMV <- paste(Yvalue, "~ 0 +", Yfunc, ":((")
  #remove any mention of timeCol
  formulaEndMV <- formulaEnd
#  formulaEndMV <- gsub(paste0(timeCol, "\\s*[:|\\*]"), "", formulaEnd)
#  formulaEndMV <- gsub(paste0("[:|\\*]\\s*", timeCol), "", formulaEndMV)
#  formulaEndMV <- gsub(paste0("\\+\\s*", timeCol), "", formulaEndMV)


  if(model == "E") # Prevent partial match of EXPR to E
  {
    list(
      fit_E(formulaStartMV, formulaEndMV, prop, ID, weightMV, funcCorr, method, dataMV, theta)[[2]],
      fit_E(formulaStartT, formulaEndT, prop, ID, weightT, timeCorr, method, dataT, theta)[[2]]
    )
  }
  else
  {
    switch(
      model,
      "STR" = list(
        fit_Intercept(MVflag = TRUE, Timeflag = FALSE, Yfunc, Yvalue, timeCol, formulaEndMV, weightMV, funcCorr, method, dataMV, exFixflag)[[2]],
        fit_Intercept(MVflag = FALSE, Timeflag = TRUE, Yfunc, Yvalue, timeCol, formulaEndT, weightT, timeCorr, method, dataT, exFixflag)[[2]]
      ),
      "ID" = list(
        fit_ID(formulaStartMV, formulaEndMV, prop, ID, weightMV, funcCorr, method, dataMV)[[2]],
        fit_ID(formulaStartT, formulaEndT, prop, ID, weightT, timeCorr, method, dataT)[[2]]
      ),
      "FULL" = list(
        fit_FULL(formulaStartMV, formulaEndMV, prop, ID, weightMV, funcCorr, method, dataMV, theta)[[2]],
        fit_FULL(formulaStartT, formulaEndT, prop, ID, weightT, timeCorr, method, dataT, theta)[[2]]
      ),
      "AV" = list(
        fit_AV(formulaStartMV, formulaEndMV, prop, ID, weightMV, funcCorr, method, dataMV, theta)[[2]],
        fit_AV(formulaStartT, formulaEndT, prop, ID, weightT, timeCorr, method, dataT, theta)[[2]]
      ),
      "ADD" = list(
        fit_ADD(formulaStartMV, formulaEndMV, prop, ID, weightMV, funcCorr, method, dataMV, theta)[[2]],
        fit_ADD(formulaStartT, formulaEndT, prop, ID, weightT, timeCorr, method, dataT, theta)[[2]]
      ),
      "FG" = list(
        fit_FG(formulaStartMV, formulaEndMV, prop, ID, weightMV, funcCorr, method, dataMV, theta)[[2]],
        fit_FG(formulaStartT, formulaEndT, prop, ID, weightT, timeCorr, method, dataT, theta)[[2]]
      ),
    )
  }
}
