
installAndLoadRequiredPackages <- function() {
  packages <- c(
    "munsell",
    "ggplot2",
    "plotly",
    "dplyr",
    "bslib",
    "shiny",
    "shinythemes",
    "httr",
    "shinylive"
  )
  
  # Function to install and load packages
  install_and_load <- function(package) {
    if (!requireNamespace(package, quietly = TRUE)) {
      install.packages(
        package,
        character.only = TRUE,
        dependencies = TRUE,
        repos = c(
          "http://cran.rstudio.com/",
          "https://cloud.r-project.org"
        )
      )
    }
    suppressMessages(library(package, character.only = TRUE))
  }
  
  
  # Install and load CRAN packages
  invisible(lapply(packages, install_and_load))
}

# installAndLoadRequiredPackages()
