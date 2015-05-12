\name{anomaly}
\alias{anomaly}
\title{Anomalous time-series detection}

\description{
  \code{anomaly} is a function for detecting unusual (i.e., anomalous, abnormal, outlier) time-series in a set of time-series.
}

\usage{
  anomaly(x, n = 10, method = c("hdr", "ahull"), robust = TRUE, plot = TRUE, 
  labels = TRUE, col)
}

\arguments{
  \item{x}{A feature matrix returned by `tsmeasures` function}
  \item{n}{Number of unusual time-series to return}
  \item{method}{Bivariate outlier detection method used for detecting high density regions of the first two principle components extracted from the time-series}
  \item{robust}{If TRUE a robust PCA will be used for feature extraction}
  \item{plot}{If TRUE, a visualization of the anomalous time-series in the first two principle compoents will be shown}
  \item{labels}{If TRUE, labels will be added to give the anomlous time series an ordering index.}
  \item{col}{A vector of length 2 giving the colours for the first and second set of points respectively (and the corresponding axes). 
             If a single colour is specified it will be used for both sets. If missing the default colour is used.}
}

\value{
  A vector of \code{n} most unusual time-series and 
  a matrix of principal component scores
} 

\author{
  Rob J Hyndman, Earo Wang, Nikolay Laptev
}

\examples{
  z <- ts(matrix(rnorm(3000),ncol = 30),freq=4)
  y <- tsmeasures(z, window = 30)
  anomaly(y, n = 2, method = "ahull") 
}
\keyword{anomaly}
