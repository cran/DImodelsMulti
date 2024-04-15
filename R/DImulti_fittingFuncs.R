###########################################################################################################################################################################################
#Internal Modelling Functions for use in DImulti_fit, DImulti_Vcov, and DImulti_theta
#' @keywords internal
#' @noRd

fit_Intercept <- function(MVflag, Timeflag, Yfunc, Yvalue, timeCol, formulaEnd, weightGen, corr, method, data, exFixflag)
{

  if(exFixflag)
  {
    formulaEndSTR <- sub(".", "", formulaEnd)
    formulaEndSTR <- sub(")[\\+:\\*]", "", formulaEndSTR)
    formulaEndSTR <- paste0(":(", formulaEndSTR)
  }
  else
  {
    formulaEndSTR <- ""
  }

  if(MVflag && Timeflag)
  {
    formula <- paste0(Yvalue, "~ 0 +", Yfunc, ":", timeCol, formulaEndSTR)
  }
  else if(MVflag)
  {
    formula <- paste0(Yvalue, "~ 0 +", Yfunc, formulaEndSTR)
  }
  else
  {
    formula <- paste0(Yvalue, "~ 0 + ", timeCol, formulaEndSTR)
  }

  formulaSTR <- stats::as.formula(formula)


  tryCatch(
    {   modelSTR <- eval(bquote( nlme::gls(.(formulaSTR), weights = .(weightGen), correlation = .(corr), method = .(method), data = data) )) },
    error = function(e) {
      message("An error has occurred")
      print(formula)
      message(e, "\n")
      stop("\nIf the error refers to the fit as 'singular' ",
           "some element of the vcov matrix is 0. Simplifying the model ",
           "can usually fix this issue") }  )

  return(list(formulaSTR, list(modelSTR)))
}


fit_ID <- function(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data)
{
  formula <- formulaStart

  #IDs
  formula <- paste0(formula, paste(ID, collapse = "+"), formulaEnd)

  formulaID <- stats::as.formula(formula)

  tryCatch(
    { modelID <- eval(bquote( nlme::gls(.(formulaID), weights = .(weightGen), correlation = .(corr), method = .(method),
                                        control = nlme::glsControl(msMaxIter = 50000, opt = "optim"), data = data) )) },
    error = function(e) {
      message("An error has occurred")
      print(formula)
      message(e, "\n")
      stop("\nIf the error refers to the fit as 'singular' ",
           "some element of the vcov matrix is 0. Simplifying the model ",
           "or changing the values of theta supplied can usually fix this issue") }  )

  return(list(formulaID, list(modelID)))
}


fit_FULL <- function(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data, theta)
{
  formula <- formulaStart

  #IDs + Pairs
  formula <- paste0(formula, paste(ID, collapse = "+"), "+", paste(colnames(dplyr::select(data, dplyr::starts_with("FULL."))), collapse = "+"),
                    formulaEnd)

  formulaPair <- stats::as.formula(formula)

  tryCatch(
    { modelPair <- eval(bquote( nlme::gls(.(formulaPair), weights = .(weightGen), correlation = .(corr), method = .(method),
                                          control = nlme::glsControl(msMaxIter = 50000, opt = "optim"), data = data) )) },
    error = function(e) {
      message("An error has occurred")
      print(formula)
      cat("Theta: ", paste(names(theta), round(theta, 4), collapse = ", "))
      message(e, "\n")
      stop("\nIf the error refers to the fit as 'singular' ",
           "some element of the vcov matrix is 0. Simplifying the model ",
           "or changing the values of theta supplied can usually fix this issue") }
  )

  return(list(formulaPair, list(modelPair)))
}


fit_E <- function(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data, theta)
{
  formula <- formulaStart

  #IDs + E
  formula <- paste0(formula, paste(ID, collapse = "+"), " + E", formulaEnd)

  formulaE <- stats::as.formula(formula)

  tryCatch(
    { modelE <- eval(bquote( nlme::gls(.(formulaE), weights = .(weightGen), correlation = .(corr), method = .(method),
                                       control = nlme::glsControl(msMaxIter = 50000, opt = "optim"), data = data) )) },
    error = function(e) {
      message("An error has occurred")
      print(formula)
      cat("Theta: ", paste(names(theta), round(theta, 4), collapse = ", "))
      message(e, "\n")
      stop("\nIf the error refers to the fit as 'singular' ",
           "some element of the vcov matrix is 0. Simplifying the model ",
           "or changing the values of theta supplied can usually fix this issue") }  )


  return(list(formulaE, list(modelE)))
}


fit_AV <- function(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data, theta)
{
  formula <- formulaStart

  #IDs + AV
  formula <- paste0(formula, paste(ID, collapse = "+"), " + AV", formulaEnd)

  formulaAV <- stats::as.formula(formula)

  tryCatch(
    { modelAV <- eval(bquote( nlme::gls(.(formulaAV), weights = .(weightGen), correlation = .(corr), method = .(method),
                                        control = nlme::glsControl(msMaxIter = 50000, opt = "optim"), data = data) )) },
    error = function(e) {
      message("An error has occurred")
      print(formula)
      cat("Theta: ", paste(names(theta), round(theta, 4), collapse = ", "))
      message(e, "\n")
      stop("\nIf the error refers to the fit as 'singular' ",
           "some element of the vcov matrix is 0. Simplifying the model ",
           "or changing the values of theta supplied can usually fix this issue") }  )

  return(list(formulaAV, list(modelAV)))
}


fit_ADD <- function(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data, theta)
{
  formula <- formulaStart

  #IDs + ADDs
  formula <- paste0(formula, paste(ID, collapse = "+"), "+", paste(colnames(dplyr::select(data, dplyr::ends_with("_add"))), collapse = "+"),
                    formulaEnd)

  formulaADD <- stats::as.formula(formula)

  tryCatch(
    { modelADD <- eval(bquote( nlme::gls(.(formulaADD), weights = .(weightGen), correlation = .(corr), method = .(method),
                                         control = nlme::glsControl(msMaxIter = 50000, opt = "optim"), data = data) )) },
    error = function(e) {
      message("An error has occurred")
      print(formula)
      cat("Theta: ", paste(names(theta), round(theta, 4), collapse = ", "))
      message(e, "\n")
      stop("\nIf the error refers to the fit as 'singular' ",
           "some element of the vcov matrix is 0. Simplifying the model ",
           "or changing the values of theta supplied can usually fix this issue") }  )

  return(list(formulaADD, list(modelADD)))
}


fit_FG <- function(formulaStart, formulaEnd, prop, ID, weightGen, corr, method, data, theta)
{
  formula <- formulaStart

  #IDs + FGs
  formula <- paste0(formula, paste(ID, collapse = "+"), "+", paste(colnames(dplyr::select(data, dplyr::starts_with("FG."))), collapse = "+"),
                    formulaEnd)

  formulaFG <- stats::as.formula(formula)

  tryCatch(
  { modelFG <- eval(bquote( nlme::gls(.(formulaFG), weights = .(weightGen), correlation = .(corr), method = .(method),
                                    control = nlme::glsControl(msMaxIter = 50000, opt = "optim"), data = data) )) },
  error = function(e) {
    message("An error has occurred")
    print(formula)
    cat("Theta: ", paste(names(theta), round(theta, 4), collapse = ", "))
    message(e, "\n")
    stop("\nIf the error refers to the fit as 'singular' ",
         "some element of the vcov matrix is 0. Simplifying the model ",
         "or changing the values of theta supplied can usually fix this issue") }  )

  return(list(formulaFG, list(modelFG)))
}
