\name{tsmeasures}
\alias{tsmeasures}
\title{Computes a feature matrix of a set of time-series}
\description{
  Computes various different features (e.g., entropy, trend, seasonality) of a set of time-series
}
\usage{
tsmeasures(y, normalise = TRUE, 
  width=ifelse(frequency(y)>1, frequency(y), 10), window=width) 
}
\arguments{
  \item{y}{A multivariate time series}
  \item{normalise}{If \code{TRUE}, each time series is scaled to be normally distributed with mean 0 and sd 1}
  \item{width}{A window size for variance change, level shift and lumpiness}
  \item{window}{A window size for KLscore}
}

\value{
  An object of class \code{features} with the following components:

  \item{lumpiness}{Variance of annual variances of remainder}
  \item{entropy}{Spectral entropy from \code{ForeCA} package}
  \item{ACF1}{First order of autocorrelation}
  \item{lshift}{Level shift using rolling window}
  \item{vchange}{Variance change}
  \item{cpoints}{The number of crossing points}
  \item{fspots}{Flat spots using disretization}
  \item{trend}{Strength of trend}
  \item{linearity}{Strength of linearity}
  \item{curvature}{Strength of curvature}
  \item{spikiness}{Strength of spikiness}
  \item{season}{Strength of seasonality}
  \item{peak}{Strength of peaks}
  \item{trough}{Strength of trough}
  \item{KLscore}{Kullback-Leibler score}
  \item{change.idx}{Index of the maximum KL score}
}

\author{
  Rob J Hyndman, Earo Wang, Nikolay Laptev
}

\seealso{
  \code{\link[ForeCA]{continuous_entropy}}
}
\examples{
  y <- tsmeasures(dat0, window = 48)
  # Add mean and variance to the "features" matrix
  avg <- colMeans(dat0)
  variance <- apply(dat0, 2, var)
  y <- cbind(y, avg, variance)
  class(y) <- c("features", "matrix")
}
\keyword{ts}
