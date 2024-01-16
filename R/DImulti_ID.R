###########################################################################################################################################################################################
#Internal ID grouping functions for use in DImulti_fit and predict.DImulti
#Code taken from DImodels package, all credit to their authors
#' @keywords internal
#' @noRd

group_IDs <- function (data, prop, ID)
{
  props <- data[, prop]
  uIDs <- unique(ID)
  nIDs <- length(uIDs)
  sapply(uIDs, function(x)
  {
    idx <- which(ID == x)
    if (length(idx) == 1)
    {
      result <- props[, idx]
    }
    else
    {
      result <- rowSums(props[, idx])
    }
  })
}


ID_name_check <- function (ID, prop, FG)
{
  cond1 <- length(grep(":", ID)) > 0
  cond2 <- any(ID == "_")
  cond3 <- any(ID == "i")
  cond4 <- any(ID == "n")
  cond5 <- any(ID == "f")
  cond6 <- any(ID == "g")
  cond7 <- any(ID == "_i")
  cond8 <- any(ID == "in")
  cond9 <- any(ID == "nf")
  cond10 <- any(ID == "fg")
  cond11 <- any(ID == "g_")
  cond12 <- any(ID == "_in")
  cond13 <- any(ID == "inf")
  cond14 <- any(ID == "nfg")
  cond15 <- any(ID == "fg_")
  cond16 <- any(ID == "_inf")
  cond17 <- any(ID == "infg")
  cond18 <- any(ID == "nfg_")
  cond19 <- any(ID == "_infg")
  cond20 <- any(ID == "infg_")
  cond21 <- any(ID == "_infg_")
  cond22 <- length(grep(pattern = "^add", x = ID, ignore.case = TRUE)) >
    0
  cond23 <- length(grep(pattern = "^full", x = ID, ignore.case = TRUE)) >
    0
  cond24 <- length(grep(pattern = "^fg", x = ID, ignore.case = TRUE)) >
    0
  cond25 <- any(ID == "E")
  cond26 <- any(ID == "AV")
  if (cond1 | cond2 | cond3 | cond4 | cond5 | cond6 | cond7 |
      cond8 | cond9 | cond10 | cond11 | cond12 | cond13 |
      cond14 | cond15 | cond16 | cond17 | cond18 | cond19 |
      cond20 | cond21 | cond22 | cond23 | cond24 | cond25 |
      cond26) {
    stop("Please give your IDs a different name.", " Names should not include colons (':'), or any single or multiple",
         " character combination of the expressions '_infg_', 'ADD', 'FG', 'FULL', 'E', and 'AV'.",
         " These expressions are reserved for internal computations.")
  }
  if (length(ID) != length(prop)) {
    stop("'ID' should be a vector of same length as prop vector specifying the grouping structure for the ID effects")
  }
  if (any(ID %in% FG)) {
    stop("'ID' cannot have any names common with names of functional groups specified by 'FG'. Change the name of groups in either 'ID' or 'FG'")
  }
}
