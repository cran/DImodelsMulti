##################################################################################################################################################################

##################################################################################################################################################################
#' predict.DImulti
#'
#' @method predict DImulti
#'
#' @description Predict from a multivariate repeated measures DI model
#'
#' @param object an object of class DImulti
#' @param newdata an optional dataframe containing the communities from which to predict. If data is
#' multivariate and in a wide format, to predict from a subset of ecosystem functions, as opposed
#' to all, please include a column for each function with any numerical value, e.g.
#' \code{newdata$Y2 <- 0}. If predicting from all functions, these columns may be included or left
#' out.
#' @param stacked a logical value used to determine whether the output is in a wide or stacked
#' format. Defaults to TRUE, meaning output is stacked/long. \cr If set to FALSE, non-unique groups of
#' unit_IDs, ecosystem function, and time points will be aggregated upon widening using the mean
#' function, please use unique unit_IDs values through newdata to avoid aggregation.
#' @param ... some methods for this generic function require additional arguments. None are used in
#' this method.
#'
#' @return The predictions from the supplied fitted DI models for the provided 'newdata', or the
#' data used to fit the model if no 'newdata' is supplied. Predictions are returned in either a
#' stacked or wide dataframe format.
#'
#' @seealso
#' \code{\link[nlme]{predict.gls}} which this function wraps.
#'
#' @export

predict.DImulti <- function(object, newdata = NULL, stacked = TRUE, ...)
{

  #Terms in model (removing response value Y)
  objTerms <- names(attr(stats::terms(object), "dataClasses"))[-1]

  if(is.null(newdata))
  {
    objData <- TRUE
    newdata <- attr(object, "data")
  }
  else if(attr(object, "MVflag"))
  {
    objData <- FALSE

    if(inherits(newdata, 'tbl_df'))
    {
      newdata <- as.data.frame(newdata)
    }

    #Missing unit_IDs column
    if(is.null(newdata[[attr(object, "unitIDs")]]))
    {
      warning("The column containing unit_IDs has not been supplied through newdata. ",
              "This column is required as a grouping factor for the covarying responses, ",
              "although its value does not matter as there is no between subject effect included. ",
              "Defaulting to row numbers. ")

      newdata[[attr(object, "unitIDs")]] <- 1:nrow(newdata)
    }

    if(attr(object, "Timeflag") & !(attr(object, "time") %in% colnames(newdata)))
    {
      newdata <- merge(newdata, unique(attr(object, "data")[[attr(object, "time")]]))
      colnames(newdata)[length(colnames(newdata))] <- attr(object, "time")
      newdata[[attr(object, "time")]] <- as.factor(newdata[[attr(object, "time")]])
    }

    #If wide data, stack
    if(length(attr(object, "y")) != 1)
    {
      if(!any(attr(object, "y") %in% colnames(newdata)))
      {
        newdata[attr(object, "y")] <- 0
      }

      predictors <- setdiff(colnames(newdata), attr(object, "y"))

      newdata <- reshape2::melt(newdata,
                                id.vars = predictors,
                                variable.name = "func",
                                value.name = "value")
    }
    else #stacked
    {
      if(!(attr(object, "Yfunc") %in% colnames(newdata)))
      {
        newdata <- merge(newdata, unique(attr(object, "data")[[attr(object, "Yfunc")]]))
        colnames(newdata)[length(colnames(newdata))] <- attr(object, "Yfunc")
      }
    }

    newdata <- newdata[order(newdata[[attr(object, "unitIDs")]], newdata[[attr(object, "Yfunc")]]),]
  }
  else # Repeated Measures Only
  {
    objData <- FALSE

    #Missing unit_IDs column
    if(is.null(newdata[[attr(object, "unitIDs")]]))
    {
      warning("The column containing unit_IDs has not been supplied through newdata. ",
              "This column is required as a grouping factor for the covarying responses, ",
              "although its value does not matter as there is no between subject effect included. ",
              "Defaulting to row numbers. ")

      newdata[[attr(object, "unitIDs")]] <- 1:nrow(newdata)
    }

    if(!(attr(object, "time") %in% colnames(newdata)))
    {
      newdata <- merge(newdata, unique(attr(object, "data")[[attr(object, "time")]]))
      colnames(newdata)[length(colnames(newdata))] <- attr(object, "time")
      newdata[[attr(object, "time")]] <- as.factor(newdata[[attr(object, "time")]])
    }

    newdata <- newdata[order(newdata[[attr(object, "unitIDs")]]),]
  }

  singleRow <- FALSE

  if(nrow(newdata) == 1)
  {
    singleRow <- TRUE
    newdata <- rbind(newdata, newdata)
  }


  missing <- c()
  #Check if all needed props are included
  if(!all(attr(object, "props") %in% colnames(newdata)))
  {
    for(i in attr(object, "props"))
    {
      if(!(i %in% colnames(newdata)))
      {
        missing <- append(missing, i)
      }
    }
    stop("The following initial proportion columns are missing from 'newdata':\n",
         paste0(missing, collapse = ", "))
  }


  if((attr(object, "DImodel") != "STR") & (attr(object, "DImodel") != "ID") & !objData)
  {

    if(inherits(newdata, 'tbl_df'))
    {
      newdata <- as.data.frame(newdata)
    }
    if(any(c("AV", "E", "NULL", "NA") %in% colnames(newdata)))
    {
      stop("You must not have any columns named \"AV\" or \"E\", \"NA\" or \"NULL\" in your dataset provided through the parameter 'newdata'")
    }
    if(any(startsWith(colnames(newdata), c("FULL.", "FG."))))
    {
      stop("You must not have any column names beginning with \"FULL.\" or \"FG.\" in your dataset provided through the parameter 'newdata'")
    }
    if(any(endsWith(colnames(newdata), c("_add"))))
    {
      stop("You must not have any column names ending with \"_add\" in your dataset provided through the parameter 'newdata'")
    }

    if(length(unique(attr(object, "thetas")) == 1)) # All the same theta value
    {
      intCols <- DImodels::DI_data(prop = attr(object, "props"), FG = attr(object, "FGs"), data = newdata, theta = attr(object, "thetas")[1],
                                   what = attr(object, "DImodel"))

      #Change new column names
      if(attr(object, "DImodel") == "FULL")
      {
        for(i in 1:ncol(intCols))
        {
          colnames(intCols)[i] <- paste0("FULL.", colnames(intCols)[i])
        }
        newdata <- cbind(newdata, data.frame(intCols))
      }
      else if(attr(object, "DImodel") == "FG")
      {
        for(i in 1:ncol(intCols))
        {
          colnames(intCols)[i] <- paste0("FG.", colnames(intCols)[i])
        }
        newdata <- cbind(newdata, data.frame(intCols))
      }
      else if(attr(object, "DImodel") %in% c("AV", "E"))
      {
        newdata <- cbind(newdata, data.frame(intCols))

        colnames(newdata)[ncol(newdata)] <- attr(object, "DImodel")
      }
      else if(attr(object, "DImodel") == "ADD")
      {
        newdata <- cbind(newdata, data.frame(intCols))
      }
    }

    else # Differing theta values
    {
      dataTemp <- data.frame()
      iCount <- 1

      #need to divide up dataset by EF and apply each theta in loop
      for(i in unique(newdata[, attr(object, "Yfunc")]))
      {
        intCols <- DImodels::DI_data(prop = attr(object, "props"), FG = attr(object, "FGs"),
                                     data = newdata[which(newdata[, attr(object, "Yfunc")] == i), ], theta = attr(object, "thetas")[iCount],
                                     what = attr(object, "DImodel"))

        #Change new column names
        if(attr(object, "DImodel") == "FULL")
        {
          for(j in 1:ncol(intCols))
          {
            colnames(intCols)[j] <- paste0("FULL.", gsub(":", ".", colnames(intCols)[j]))
          }
        }
        else if(attr(object, "DImodel") == "FG")
        {
          for(j in 1:ncol(intCols))
          {
            colnames(intCols)[j] <- paste0("FG.", colnames(intCols)[j])
          }
        }

        dataTemp <- rbind(dataTemp, cbind(newdata[which(newdata[, attr(object, "Yfunc")] == i), ], intCols))

        iCount <- iCount + 1
      }
      newdata <- dataTemp

      if(attr(object, "DImodel") %in% c("AV", "E"))
      {
        colnames(newdata)[ncol(newdata)] <- attr(object, "DImodel")
      }

      newdata <- newdata[order(newdata[[attr(object, "unitIDs")]], newdata[[attr(object, "Yfunc")]]), ]
    }
  }

  #Non unique rows which would be aggregated upon widening
  # if(!all(!duplicated(newdata[, !(names(newdata) %in% attr(object, "Yvalue"))])) & #!stacked &
  #    !singleRow)
  # {
  #   newdata <- dplyr::distinct(newdata[, !(names(newdata) %in% attr(object, "Yvalue"))])
  #   warning("Removing rows with duplicate predictors in newdata, ",
  #           "you can use function duplicated() to test your dataset")
  #
  #   print(newdata)
  # }



  #ID grouping
  ID_name_check(ID = attr(object, "IDs"), prop = attr(object, "props"), FG = attr(object, "FGs"))
  grouped_ID <- group_IDs(data = newdata, prop = attr(object, "props"), ID = attr(object, "IDs"))

  newdata <- cbind(newdata, grouped_ID)

  missing <- c()
  #Check if all needed columns are included
  if(!all(objTerms %in% colnames(newdata)))
  {
    for(i in objTerms)
    {
      if(!(i %in% colnames(newdata)))
      {
        if(is.factor(attr(object, "data")[[i]]))
        {
          warning("The model term ", i, " is missing from  the dataset supplied through newdata. ",
                  "The column has been added with a set value taken from the training dataset: ",
                  levels(attr(object, "data")[[i]])[1])
          newdata[[i]] <- levels(attr(object, "data")[[i]])[1]
        }
        else if(is.numeric(attr(object, "data")[[i]]))
        {
          warning("The model term ", i, " is missing from  the dataset supplied through newdata. ",
                  "The column has been added with a mean value taken from the training dataset: ",
                  mean(attr(object, "data")[[i]]))
          newdata[[i]] <- mean(attr(object, "data")[[i]])
        }
        else #fail safe for non-numeric non-factor columns
        {
          missing <- append(missing, i)
        }
      }
    }

    if(!is.null(missing)) # missing non-factor column
    {
      stop("The following model terms are missing from the dataset supplied through newdata:\n",
           paste0(i, collapse = ", "))
    }
  }


  #Number of responses to be predicted
  nFuncs <- 1
  nTimes <- 1

  if(attr(object, "MVflag"))
  {
    nFuncs <- nlevels(as.factor(newdata[, attr(object, "Yfunc")]))
  }
  if(attr(object, "Timeflag"))
  {
    nTimes <- nlevels(as.factor(newdata[, attr(object, "time")]))
  }

  #Predict for each community, add each row to modelPreds (dataframe)
  form <- stats::delete.response(object[["terms"]])
  contr <- object$contrasts
  dataMod <- stats::model.frame(formula = form, data = newdata,
                                drop.unused.levels = TRUE, xlev = lapply(contr, rownames))
  N <- nrow(dataMod)
  if (length(all.vars(form)) > 0)
  {
    X <- stats::model.matrix(form, dataMod, contr)
  }
  else
  {
    X <- array(1, c(N, 1), list(row.names(dataMod), "(Intercept)"))
  }
  cf <- object$coefficients
  val <- c(X[, names(cf), drop = FALSE] %*% cf)
  lab <- "Predicted values"
  if (!is.null(aux <- attr(object, "units")$y))
  {
    lab <- paste(lab, aux)
  }
  modelPreds <- structure(val, label = lab)

  if(singleRow)
  {
    modelPreds <- modelPreds[1]
  }

  if(attr(object, "MVflag") & !attr(object, "Timeflag"))
  {
    modelPreds <- cbind.data.frame(newdata[[attr(object, "unitIDs")]], modelPreds, newdata[, attr(object, "Yfunc")])
    colnames(modelPreds) <- c(attr(object, "unitIDs"), "Yvalue", attr(object, "Yfunc"))
  }
  else if(attr(object, "Timeflag") & !attr(object, "MVflag"))
  {
    modelPreds <- cbind.data.frame(newdata[[attr(object, "unitIDs")]], modelPreds, newdata[, attr(object, "time")])
    colnames(modelPreds) <- c(attr(object, "unitIDs"), "Yvalue", attr(object, "time"))
  }
  else #both
  {
    modelPreds <- cbind.data.frame(newdata[[attr(object, "unitIDs")]], modelPreds)
    modelPreds$Ytype <- apply(newdata[, c(attr(object, "Yfunc"), attr(object, "time"))], 1, paste, collapse = ":")
    colnames(modelPreds) <- c(attr(object, "unitIDs"), "Yvalue", "Ytype")
  }

  if(singleRow)
  {
    modelPreds <- modelPreds[1, ]
  }

  if(!stacked)
  {
    predForm <- stats::as.formula(paste(attr(object, "unitIDs"), "~", colnames(modelPreds)[3]))

    modelPreds <- reshape2::dcast(data = modelPreds, formula = predForm, value.var = "Yvalue", fun.aggregate = mean)
  }

  return(modelPreds)
}
