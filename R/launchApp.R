#' launchApp
#'
#' Run OncoCast App, enabling the user to either generate an OncoCast run (though not recommended)
#' if it can be ran in the console. Or create a flexible framework to explore the results of an
#' OncoCast run that has been saved.
#' @export
#' @import
#' shiny
#' shinydashboard

launchApp <- function(){

  runApp(
    appDir = find.package("OncoCast")
    # appDir = paste0("R")
  )
}

