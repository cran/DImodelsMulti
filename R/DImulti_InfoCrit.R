##################################################################################################################################################################
#' BICc
#'
#' @description Calculate Second-order Bayesian Information Criterion for a fitted DImulti obect
#'
#' @param model an object of class DImulti
#'
#' @return A numeric value of the BICc value corresponding to the model object
#'
#' @export

BICc <- function(model)
{
  BIC <- BIC(model)
  n <- length(model$fitted)
  K <- attr(stats::logLik(model), "df")

  BICc <- BIC + (log(n) * (K + 1) * K)/(n - K - 1)

  BICc
}

##################################################################################################################################################################

##################################################################################################################################################################
#' AICc
#'
#' @description Calculate Second-order Akaike Information Criterion for a fitted DImulti object
#'
#' @param model an object of class DImulti
#'
#' @return A numeric value of the AICc value corresponding to the model object
#'
#' @export

AICc <- function(model)
{
  AIC <- AIC(model)
  n <- length(model$fitted)
  K <- attr(stats::logLik(model), "df")

  AICc <- AIC + (2 * K^2 + 2 * n)/(n - K - 1)

  AICc
}

##################################################################################################################################################################
