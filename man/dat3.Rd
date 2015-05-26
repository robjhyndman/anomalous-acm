\name{dat3}
\alias{dat3}
\title{Aggregated and anonymized datasets from Yahoo representing server metrics of Yahoo services}
\description{
  Aggregated and anonymized dataset from Yahoo representing server metrics of Yahoo services.
  \code{dat3} contains 571 time series each with 741 observations.
}

\examples{
  plot(dat3[,1:10])
  plot(dat3, plot.type='single', col=1:50)
}

\format{time series}

\author{
  Rob J Hyndman, Earo Wang, Nikolay Laptev
}

\keyword{datasets}
