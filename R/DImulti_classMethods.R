##################################################################################################################################################################
#' print.DImulti
#' @method print DImulti
#'
#' @description Print details of the fitted DI models supplied
#'
#' @param x an object of class DImulti
#' @param ... some methods for this generic function require additional arguments. None are used in
#' this method.
#'
#' @return object x
#'
#' @details
#' The appearance of the printed information will differ depending on whether a user has installed
#' some combination of the suggested packages "crayon", "cli", and "fansi".
#' These changes are mainly cosmetic, with crayon making the output easier to read, cli providing
#' links to help files, and fansi enabling the reading of special characters in R markdown (Rmd)
#' files. See 'Examples' below for suggested code to include in Rmd files.
#'
#'
#' @seealso
#' \code{\link[base]{print}} which this function wraps.
#'
#'
#' @examples
#' ## Set up for R markdown for crayon and cli output if user has packages installed
#'
#' if(requireNamespace("fansi", quietly = TRUE) &
#'    requireNamespace("crayon", quietly = TRUE) &
#'    requireNamespace("cli", quietly = TRUE))
#' {
#'  options(crayon.enabled = TRUE)
#'
#'  ansi_aware_handler <- function(x, options)
#'  {
#'   paste0(
#'     "<pre class=\"r-output\"><code>",
#'     fansi::sgr_to_html(x = x, warn = FALSE, term.cap = "256"),
#'     "</code></pre>"
#'   )
#'  }
#'
#'  old_hooks <- fansi::set_knit_hooks(knitr::knit_hooks,
#'                                     which = c("output", "message", "error", "warning"))
#'
#'  knitr::knit_hooks$set(
#'    output = ansi_aware_handler,
#'    message = ansi_aware_handler,
#'    warning = ansi_aware_handler,
#'    error = ansi_aware_handler
#'   )
#'  }
#' #################################################################################################
#'
#' ## Usage
#'
#' model <- DImulti(prop = c("G1", "G2", "L1", "L2"), y = "Y", eco_func = c("Var", "un"),
#'                  unit_IDs = "Plot", theta = c(0.5, 1, 1.2), DImodel = "FG",
#'                  FG = c("Grass", "Grass", "Legume", "Legume"), extra_fixed = ~ Density,
#'                  method = "REML", data = dataBEL)
#'
#' print(model)
#'
#' @export

print.DImulti <- function(x, ...)
{
  if(requireNamespace("crayon", quietly = TRUE))
  {
    #Print notes
    cat(crayon::bold("Note: \n"))
    cat(crayon::bold(paste(crayon::magenta("Method Used ="), crayon::underline(attr(x, "method")), "\n")))
    cat(crayon::bold(crayon::magenta("Correlation Structure Used = ")))

    if(cli::ansi_has_hyperlink_support() & requireNamespace("cli", quietly = TRUE))
    {
      switch(tolower(attr(x, "correlation")),
             "un" = cat(crayon::bold(crayon::underline(cli::cli_text("{.help [UN](nlme::corSymm)}")))),
             "cs" = cat(crayon::bold(crayon::underline(cli::cli_text("{.help [CS](nlme::corCompSymm)}")))),
             "ar1" = cat(crayon::bold(crayon::underline(cli::cli_text("{.help [AR(1)](nlme::corAR1)}")))),
             "un@un" = cat(crayon::bold(crayon::underline(cli::cli_text("{.help [UN](nlme::corSymm)} @ {.help [UN](nlme::corSymm)}")))),
             "un@cs" = cat(crayon::bold(crayon::underline(cli::cli_text("{.help [UN](nlme::corSymm)} @ {.help [CS](nlme::corCompSymm)}")))),
             "un@ar1" = cat(crayon::bold(crayon::underline(cli::cli_text("{.help [UN](nlme::corSymm)} @ {.help [AR(1)](nlme::corAR1)}")))),
             "cs@un" = cat(crayon::bold(crayon::underline(cli::cli_text("{.help [CS](nlme::corCompSymm)} @ {.help [UN](nlme::corSymm)}")))),
             "cs@cs" = cat(crayon::bold(crayon::underline(cli::cli_text("{.help [CS](nlme::corCompSymm)} @ {.help [CS](nlme::corCompSymm)}")))),
             "cs@ar1" = cat(crayon::bold(crayon::underline(cli::cli_text("{.help [CS](nlme::corCompSymm)} @ {.help [AR(1)](nlme::corAR1)}"))))
      )
    }
    else
    {
      cat(crayon::bold(crayon::underline(attr(x, "correlation"))))
    }

    #Print Model Name
    cat("\n", crayon::bold(crayon::green(attr(x, "name"))))
    #If an interaction model, print theta
    if(attr(x, "name") != "Intercept Only Model" & attr(x, "name") != "ID Model")
    {
      if(attr(x, "estThetas") == TRUE)
      {
        cat(crayon::bold(crayon::magenta("\nTheta estimate(s) = ")))
      }
      else
      {
        cat(crayon::bold(crayon::magenta("\nTheta value(s) = ")))
      }
      if(!is.null(names(x$theta)))
      {
        cat(paste(crayon::bold(names(x$theta)), round(x$theta, 4), sep = ":", collapse = ", "))
      }
      else
      {
        cat(round(x$theta, 4), sep = ",")
      }
    }
  }
  else
  {
    #Print notes
    cat("Note: \n")
    cat(paste("Method Used =", attr(x, "method")), "\n")
    cat("Correlation Structure Used = ")

    if(cli::ansi_has_hyperlink_support() & requireNamespace("cli", quietly = TRUE))
    {
      switch(tolower(attr(x, "correlation")),
             "un" = cat(cli::cli_text("{.help [UN](nlme::corSymm)}")),
             "cs" = cat(cli::cli_text("{.help [CS](nlme::corCompSymm)}")),
             "ar1" = cat(cli::cli_text("{.help [AR(1)](nlme::corAR1)}")),
             "un@un" = cat(cli::cli_text("{.help [UN](nlme::corSymm)} @ {.help [UN](nlme::corSymm)}")),
             "un@cs" = cat(cli::cli_text("{.help [UN](nlme::corSymm)} @ {.help [CS](nlme::corCompSymm)}")),
             "un@ar1" = cat(cli::cli_text("{.help [UN](nlme::corSymm)} @ {.help [AR(1)](nlme::corAR1)}")),
             "cs@un" = cat(cli::cli_text("{.help [CS](nlme::corCompSymm)} @ {.help [UN](nlme::corSymm)}")),
             "cs@cs" = cat(cli::cli_text("{.help [CS](nlme::corCompSymm)} @ {.help [CS](nlme::corCompSymm)}")),
             "cs@ar1" = cat(cli::cli_text("{.help [CS](nlme::corCompSymm)} @ {.help [AR(1)](nlme::corAR1)}"))
      )
    }
    else
    {
      cat(attr(x, "correlation"))
    }

    #Print Model Name
    cat("\n", attr(x, "name"))
    #If an interaction model, print theta
    if(attr(x, "name") != "Intercept Only Model" & attr(x, "name") != "ID Model")
    {
      if(attr(x, "estThetas") == TRUE)
      {
        cat("\nTheta estimate(s) = ")
      }
      else
      {
        cat("\nTheta value(s) = ")
      }
      if(!is.null(names(x$theta)))
      {
        cat(paste(names(x$theta), round(x$theta, 4), sep = ":", collapse = ", "))
      }
      else
      {
        cat(round(x$theta, 4), sep = ",")
      }
    }
  }

  cat("\n")

  dd <- x$dims
  mCall <- x$call

  if(inherits(x, "gnls"))
  {
    cat("\nGeneralized nonlinear least squares fit\n")
  }
  else
  {
    cat("\nGeneralized least squares fit by ")
    cat(if(x$method == "REML") "REML\n"
        else "maximum likelihood\n")
  }
  cat("  Model:", deparse(x$call$model), "\n")
  if(!is.null(x$call$subset))
  {
    cat("  Subset:", deparse(stats::asOneSidedFormula(x$call$subset)[[2]]),
        "\n")
  }
  print(stats::setNames(c(stats::AIC(x), stats::BIC(x), stats::logLik(x)), c("AIC", "BIC", "logLik")))

  if(attr(x, "MVflag"))
  {
    cat("\n Multivariate ")
    print(summary(x$corr$`Multivariate`))
  }
  if(attr(x, "Timeflag"))
  {
    cat("\n Repeated Measure ")
    print(summary(x$corr$`Repeated Measure`))
  }

  tTable <- invisible(summary(x, notes = FALSE)$tTable)
  stars <- ifelse(tTable[, 4] < 0.001, "***",
                  ifelse(tTable[, 4] <= 0.01, "**",
                         ifelse(tTable[, 4] <= 0.05, "*",
                                ifelse(tTable[, 4] <= 0.1, "+", "")))
  )
  if(requireNamespace("crayon", quietly = TRUE))
  {
    base::print(knitr::kable(cbind(formatC(tTable[, c(1)], format = "f", digits = 3, flag = "+"),
                                   formatC(tTable[, c(2)], format = "f", digits = 3),
                                   formatC(tTable[, c(3)], format = "f", digits = 3),
                                   formatC(tTable[, c(4)], format = "g"),
                                   stars),
                             caption = crayon::bold("Fixed Effect Coefficients"), col.names = c("Beta", "Std. Error", "t-value", "p-value", "Signif"),
                             format = "simple"))
  }
  else
  {
    base::print(knitr::kable(cbind(formatC(tTable[, c(1)], format = "f", digits = 3, flag = "+"),
                                   formatC(tTable[, c(2)], format = "f", digits = 3),
                                   formatC(tTable[, c(3)], format = "f", digits = 3),
                                   formatC(tTable[, c(4)], format = "g"),
                                   stars),
                             caption = "Fixed Effect Coefficients", col.names = c("Beta", "Std. Error", "t-value", "p-value", "Signif"),
                             format = "simple"))
  }

  cat("\n")
  cat(
    noquote(
      "Signif codes: 0-0.001 '***', 0.001-0.01 '**', 0.01-0.05 '*', 0.05-0.1 '+', 0.1-1.0 ' '"))

  cat("\n\n")

  cat("Degrees of freedom:", dd[["N"]], "total;", dd[["N"]] -
        dd[["p"]], "residual\n")
  cat("Residual standard error:", format(x$sigma), "\n\n")

  print(getVarCov(x))

  invisible(x)
}


##################################################################################################################################################################
#' summary.DImulti
#' @method summary DImulti
#'
#' @description Print a summary of the fitted DI models supplied
#'
#' @param object an object inheriting from class DImulti, representing a generalized least squared fitted linear model using the
#' Diversity-Interactions framework
#' @param verbose an optional logical value used to control the amount of output when the object is printed. Defaults to FALSE.
#' @param ... some methods for this generic function require additional arguments. None are used in
#' this method.
#'
#' @return an object inheriting from class summary.gls with all components included in object (see glsObject for a full description of the
#' components) plus the following components:
#'
#' corBeta,
#' approximate correlation matrix for the coefficients estimates
#'
#' tTable,
#' a matrix with columns Value, Std. Error, t-value, and p-value representing respectively the coefficients estimates, their approximate
#' standard errors, the ratios between the estimates and their standard errors, and the associated p-value under a
#' t approximation. Rows correspond to the different coefficients.
#'
#' residuals,
#' if more than five observations are used in the gls fit, a vector with the minimum, first quartile, median, third quartile, and maximum of
#' the residuals distribution; else the residuals.
#'
#' AIC,
#' the Akaike Information Criterion corresponding to object.
#'
#' BIC,
#' the Bayesian Information Criterion corresponding to object.
#'
#' theta,
#' the values of the non-linear parameter theta used in the model
#'
#' @seealso
#' \code{\link[nlme]{summary.gls}} which this function wraps.
#'
#' @export
#'

summary.DImulti <- function(object, verbose = FALSE, ...)
{
  fixSig <- attr(object[["modelStruct"]], "fixedSigma")
  fixSig <- !is.null(fixSig) && fixSig

  stdBeta <- sqrt(diag(as.matrix(object$varBeta)))
  corBeta <- t(object$varBeta/stdBeta)/stdBeta
  beta <- object$coefficients

  dims <- object$dims
  dimnames(corBeta) <- list(names(beta), names(beta))
  object$corBeta <- corBeta

  tTable <- data.frame(beta, stdBeta, beta/stdBeta, beta)
  dimnames(tTable) <- list(names(beta), c("Value", "Std.Error",
                                          "t-value", "p-value"))
  tTable[, "p-value"] <- 2 * stats::pt(-abs(tTable[, "t-value"]),
                                dims$N - dims$p)

  object$tTable <- as.matrix(tTable) #Save tTable to object

  resd <- stats::resid(object, type = "pearson")
  if (length(resd) > 5)
  {
    resd <- stats::quantile(resd, na.rm = TRUE)
    names(resd) <- c("Min", "Q1", "Med", "Q3", "Max")
  }
  object$residuals <- resd
  aux <- stats::logLik(object)

  structure(c(object, list(BIC = stats::BIC(aux), AIC = stats::AIC(aux))),
            verbose = verbose, class = c("summary.DImulti", "summary.gls", class(object)))
}

##################################################################################################################################################################

##################################################################################################################################################################
#' print.summary.DImulti
#'
#' @method print summary.DImulti
#'
#' @description print the details of an object of class summary.DImulti
#'
#' @param x an object of class DImulti
#' @param verbose an optional logical value used to control the amount of output when the object is printed. Defaults to FALSE.
#' @param digits the number of significant digits to use when printing
#' @param ... some methods for this generic function require additional arguments. None are used in
#' this method.
#'
#' @return object x
#'
#' @seealso
#' \code{\link[base]{print}} which this function wraps.
#'
#' @export
#'

print.summary.DImulti <- function(x, verbose = FALSE, digits = .Options$digits, ...)
{
  dd <- x$dims
  fixSig <- attr(x[["modelStruct"]], "fixedSigma")
  fixSig <- !is.null(fixSig) && fixSig
  verbose <- verbose || attr(x, "verbose")

  mCall <- x$call
  if (inherits(x, "gnls")) {
    cat("Generalized nonlinear least squares fit\n")
  }
  else {
    cat("Generalized least squares fit by ")
    cat(if (x$method == "REML")
      "REML\n"
      else "maximum likelihood\n")
  }

  cat("  Model:", deparse(mCall$model), "\n")
  cat("  Data:", deparse(mCall$data), "\n")
  if (!is.null(mCall$subset))
  {
    cat("  Subset:", deparse(stats::asOneSidedFormula(mCall$subset)[[2]]),
        "\n")
  }
  print(data.frame(AIC = x$AIC, BIC = x$BIC, logLik = as.vector(x$logLik),
                   row.names = " "), ...)
  if (verbose)
  {
    cat("Convergence at iteration:", x$numIter, "\n")
  }
  if (length(x$modelStruct))
  {
    cat("\n")
    print(summary(x$modelStruct), ...)
  }
  cat("\nCoefficients:\n")
  xtTab <- as.data.frame(x$tTable)
  wchPval <- match("p-value", names(xtTab))
  for (i in names(xtTab)[-wchPval]) {
    xtTab[, i] <- format(zapsmall(xtTab[, i]))
  }
  xtTab[, wchPval] <- format(round(xtTab[, wchPval], 4L))
  if (any(wchLv <- (as.double(levels(xtTab[, wchPval])) == 0)))
  {
    levels(xtTab[, wchPval])[wchLv] <- "<.0001"
  }
  row.names(xtTab) <- dimnames(x$tTable)[[1]]
  print(xtTab, ...)

  cat("\nTheta values: ")

  if(!is.null(names(x$theta)))
  {
    cat(paste(names(x$theta), round(x$theta, 4), sep = ":", collapse = ", "))
  }
  else
  {
    cat(round(x$theta, 4), sep = ",")
  }

  cat("\n\n")

  if (nrow(x$tTable) > 1L) {
    corr <- x$corBeta
    class(corr) <- "correlation"
    print(corr, title = "\n Correlation:", ...)
  }
  cat("\nStandardized residuals:\n")
  print(x$residuals, ...)
  cat("\n")
  cat("Residual standard error:", format(x$sigma), "\n")
  cat("Degrees of freedom:", dd[["N"]], "total;", dd[["N"]] -
        dd[["p"]], "residual\n")
  invisible(x)
}


##################################################################################################################################################################

##################################################################################################################################################################
#' getVarCov.DImulti
#'
#' @method getVarCov DImulti
#'
#' @description Extract the variance-covariance matrix from a fitted DImulti model
#'
#' @param obj an object of class DImulti
#' @param ... some methods for this generic function require additional arguments. None are used in
#' this method.
#'
#' @return A variance-covariance matrix or a named list of variance-covariance matrices.
#'
#' @seealso
#' \code{\link[DImodelsMulti]{DImulti}}
#'
#' @export
#' @importFrom nlme getVarCov
#'

getVarCov.DImulti <- function(obj, ...)
{
  if(attr(obj, "MVflag") & attr(obj, "Timeflag")) #Both
  {
    return(list("Multivariate" = obj$vcov$`Multivariate`,
                "Repeated Measure" = obj$vcov$`Repeated Measure`,
                "Combined" = obj$vcov$`Combined`))
  }

  else #Tflag or MVflag
  {
    return(obj$vcov[[1]])
  }

  return(obj$vcov[[1]]) #default case
}


