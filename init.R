packages <- c("munsell","ggplot2", "plotly", "dplyr", "bslib", "shiny", "shinythemes","shinylive","httr")

# Function to install and load packages
install_and_load <- function(package) {
  if (!requireNamespace(package, quietly = TRUE)) {
    install.packages(package, character.only = TRUE, repos ="http://cran.rstudio.com/" ) #"https://cloud.r-project.org"
  }
  suppressMessages(library(package, character.only = TRUE))
}


# Install and load CRAN packages
invisible(lapply(packages, install_and_load))