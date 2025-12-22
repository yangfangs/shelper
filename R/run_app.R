#' Run the Shelper Shiny App
#'
#' @param port The port to run the application on. Default is 3838.
#' @param launch.browser Whether to launch the browser. Default is TRUE.
#' @export
run_shelper <- function(port = 3838, launch.browser = TRUE) {
  app_dir <- system.file("shiny", "shelper_app", package = "shelper")
  if (app_dir == "") {
    stop("Could not find example directory. Try re-installing `shelper`.", call. = FALSE)
  }
  
  # Load configuration if available
  # You might want to allow passing config parameters here as well
  
  shiny::runApp(app_dir, port = port, launch.browser = launch.browser, display.mode = "normal")
}
