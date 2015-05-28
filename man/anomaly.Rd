\name{anomaly}
\alias{anomaly}
\title{Anomalous time-series detection}

\description{
  \code{anomaly} is a function for detecting unusual (i.e., anomalous, abnormal, outlier) time-series in a set of time-series.
}

\usage{
  anomaly(x, n = 10, method = c("hdr", "ahull"), robust = TRUE, ordered = FALSE,
  plot = TRUE, labels = TRUE, col)
}

\arguments{
  \item{x}{A feature matrix returned by `tsmeasures` function}
  \item{n}{Number of unusual time-series to return}
  \item{method}{Bivariate outlier detection method used for detecting high density regions of the first two principal components extracted from the time-series}
  \item{robust}{If \code{TRUE} a robust PCA will be used for feature extraction}
  \item{ordered}{If \code{TRUE}, the outliers are returned in order. This only applies when
    \code{method="ahull"} as the outliers are always ordered when \code{method="hdr"}.}
  \item{plot}{If \code{TRUE}, a visualization of the anomalous time-series in the first two principal compoents will be shown}
  \item{labels}{If \code{TRUE}, labels will be added to give the anomlous time series an ordering index.}
  \item{col}{A vector of length 2 giving the colours for the first and second set of points respectively (and the corresponding axes). 
             If a single colour is specified it will be used for both sets. If missing the default colour is used.}
}

\value{
  A vector showing the row index of the \code{n} most unusual time-series and 
  a matrix of principal component scores. Returned invisibly if \code{plot=TRUE}.
} 

\details{When method = "hdr", the outliers are those with the lowest density in the 2-dimensional PC space using a bivariate kernel density estimator. In this case, the outliers are always in order of increasing density. 

When method = "ahull", the outliers are the first \code{n} scores in the 2-dimensional PC space which break away as singletons when alpha-convex hulls are used with decreasing alpha. They are not ordered by default. If \code{ordered=TRUE}, then the outliers are ordered based on when they first break away from the remaining points. This is much slower, especially for large \code{n}, and so it is turned off by default.
}

\author{
  Rob J Hyndman, Earo Wang, Nikolay Laptev
}

\seealso{
  \code{\link[alphahull]{ahull}}, \code{\link[hdrcde]{hdr.2d}}, 
}
\examples{
  y <- tsmeasures(dat0, window = 24)
  anomaly(y, n = 2, method = "ahull", ordered=TRUE) 
}
\keyword{anomaly}
