\name{biplot}
\alias{biplot.features}
\title{\code{biplot} of (robust) PCA components of the feature matrix}
\description{Plot a biplot on the current graphics device of the (robust) PCA components generated from the features returned by \code{tsmeasures}}

\arguments{
  \item{x}{Feature matrix returned by \code{tsmeasures}}
  \item{robust}{If TRUE, a robust PCA will be used on the feature matrix}
  \item{col}{A vector of length 2 giving the colours for the first and second set of points respectively (and the corresponding axes). 
             If a single colour is specified it will be used for both sets. If missing the default colour is used.}
  \item{...}{Other standard options of \code{biplot}}
}

\value{
  a plot is produced on the current graphics device.
}

\author{
  Rob J Hyndman, Earo Wang, Nikolay Laptev
}

\examples{
  z <- ts(matrix(rnorm(3000),ncol=100),freq=4)
  y <- tsmeasures(z)
  biplot(y)
}
\keyword{ts}
