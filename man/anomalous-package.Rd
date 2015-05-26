\name{anomalous-package}
\alias{anomalous-package}
\alias{anomalous}
\docType{package}
\title{Unusual Time Series Detection}

\description{
  It is becoming increasingly common for organizations to collect very large amounts of data over time, and to need to detect unusual or anomalous time series.
  For example, Yahoo has banks of mail servers that are monitored over time. Many measurements on server performance are collected every hour for each of thousands of servers. 
  A common use-case is to identify servers that are behaving unusually.
  Methods in this package compute a vector of features on each time series, measuring characteristics of the series. For example, the features may include lag correlation, strength of seasonality, 
  spectral entropy, etc. Then a principal component decomposition is used on the features, and various bivariate outlier detection methods are applied to the first two principal components.
  This enables the most unusual series, based on their feature vectors, to be identified. The bivariate outlier detection methods used are based on highest density regions and alpha-hulls.
  This package contains both synthetic and real data from Yahoo.
}
\author{
  Rob J Hyndman, Earo Wang, Nikolay Laptev
  Maintainer: Earo Wang <earo.wang@gmail.com>
}

\keyword{
  time series
  anomaly detection
  unsual time series detection
  anomalous time series detection
}

\examples{
  y <- tsmeasures(dat0)
  biplot(y)
  anomaly(y)
  anomaly(y, method="ahull")
}
