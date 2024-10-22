.onLoad <- function(libname, pkgname) {
  suppressPackageStartupMessages({
    library(fastDummies, quietly = TRUE)
    # Repeat for other libraries, even if they are in NAMESPACE
  })
}
