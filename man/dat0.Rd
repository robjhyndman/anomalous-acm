\name{dat0}
\alias{dat0}
\title{Aggregated and anonymized datasets from Yahoo representing server metrics of Yahoo services}
\description{
  Aggregated and anonymized dataset from Yahoo representing server metrics of Yahoo services.
  \code{dat0} contains 44 time series each with 1437 observations.
}

\examples{
  plot(dat0[,1:10])
  plot(dat0, plot.type='single', col=1:44)
}

\format{time series}

\author{
  Rob J Hyndman, Earo Wang, Nikolay Laptev
}

\keyword{datasets}
