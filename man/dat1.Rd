\name{dat1}
\alias{dat1}
\title{Aggregated and anonymized datasets from Yahoo representing server metrics of Yahoo services}
\description{
  Aggregated and anonymized dataset from Yahoo representing server metrics of Yahoo services.
  \code{dat1} contains 571 time series each with 741 observations.
}

\examples{
  plot(dat1[,11:20])
  plot(dat1, plot.type='single', col=1:50)
}

\format{time series}

\author{
  Rob J Hyndman, Earo Wang, Nikolay Laptev
}

\keyword{datasets}
