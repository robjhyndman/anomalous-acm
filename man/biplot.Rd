\name{biplot}
\alias{biplot.features}
\title{\code{biplot} of (robust) PCA components of the feature matrix}
\description{Plot a biplot on the current graphics device of the (robust) PCA components generated from the features returned by \code{tsmeasures}}

\usage{\method{biplot}{features}(x, robust = TRUE, col, ...)}
\arguments{
  \item{x}{Feature matrix returned by \code{tsmeasures}}
  \item{robust}{If TRUE, a robust PCA will be used on the feature matrix}
  \item{col}{A vector of length 2 giving the colours for the first and second set of points respectively (and the corresponding axes). 
             If a single colour is specified it will be used for both sets. If missing the default colour is used.}
  \item{...}{Other standard options of \code{\link[stats]{biplot}}.}
}

\value{
  a plot is produced on the current graphics device.
}

\seealso{
  \code{\link[stats]{biplot}}, \code{\link{tsmeasures}}.
}

\author{
  Rob J Hyndman, Earo Wang, Nikolay Laptev
}

\examples{
  y <- tsmeasures(dat0)
  biplot(y)
}
\keyword{ts}
