.onAttach <- function(...) {
  version <- library(help=anomalousACM)$info[[1]]
  version <- version[pmatch("Version", version)]
  um <- strsplit(version, " ")[[1]]
  version <- um[nchar(um) > 0][2]
  packageStartupMessage(paste("This is anomalousACM", version, "\n"))
}
