.onAttach <- function(...) {
  version <- library(help=anomalous)$info[[1]]
  version <- version[pmatch("Version", version)]
  um <- strsplit(version, " ")[[1]]
  version <- um[nchar(um) > 0][2]
  packageStartupMessage(paste("This is anomalous", version, "\n"))
}
