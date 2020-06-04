# SPECIFY PACKAGES -------------------------------------------------------------

packages <- c(
  "plyr", 
  "dplyr", 
  "tidyr", 
  "stringr", 
  "rlang", 
  "readr", 
  "parallel", 
  "getopt", 
  "readxl")

# FUNCTIONS --------------------------------------------------------------------

pkg_test <- function(x) {
  if ( ! x %in% rownames(installed.packages()) ) {
    # if a package isn't installed, notify the user and change exit code to 1
    message(paste0('\nERROR:\tR package "', x, '" not installed'))
    exit_code <<- 1
  }
}

# MAIN ------------------------------------------------------------------------

exit_code <- 0

# check that packages are installed
for ( package in packages ) {
  pkg_test(package)
}

# if any of the packages are not installed, exit with status 1
if ( ! exit_code == 0 ) {
  quit(save = "no", status = 1)
} else {
  quit(save = "no", status = 0)
}
