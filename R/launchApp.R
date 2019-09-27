#' launchApp
#' @export

launchApp <- function() {
  shinyApp(ui = ui, server = server)
}
