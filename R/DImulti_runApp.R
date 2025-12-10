##################################################################################################################################################################
#' DImultiApp
#'
#' @description A function to launch the DImulti R Shiny App from within an R
#' session.
#'
#'
#' @export
DImultiApp <- function() {
  # Packages needed for shiny app
  pkgs <- c("ggplot2", "DImodelsVis", "shinydashboard", "bslib")
  # Packages to be installed
  need_install <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

  if (length(need_install) > 0) {
    stop(
      paste0(
        "The following package",
        if (length(need_install) > 1) "s",
        " must be installed to run the DImodelsMulti app: ",
        paste(shQuote(need_install), collapse = ", "),"\n",
        "Install", if (length(need_install) > 1) " them" else " it",
        " with: install.packages(c(",
        paste(shQuote(need_install), collapse = ", "), "))"
      ),
      call. = FALSE
    )
  }

  #shiny::runApp("inst/shiny/DImultiApp.R")
  shiny::runApp(file.path(system.file(package = "DImodelsMulti"), "shiny", "DImultiApp.R"))
}
