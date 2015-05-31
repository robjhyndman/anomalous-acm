Anomalous time-series R Package
===================

It is becoming increasingly common for organizations to collect very large amounts of data over time, and to need to detect unusual or anomalous time series.
For example, Yahoo has banks of mail servers that are monitored over time. Many measurements on server performance are collected every hour for each of thousands of servers. 
A common use-case is to identify servers that are behaving unusually.
Methods in this package compute a vector of features on each time series, measuring characteristics of the series. For example, the features may include lag correlation, strength of seasonality, 
spectral entropy, etc. Then a robust principal component decomposition is used on the features, and various bivariate outlier detection methods are applied to the first two principal components.
This enables the most unusual series, based on their feature vectors, to be identified. The bivariate outlier detection methods used are based on highest density regions and alpha-hulls.
For demo purposes, this package contains both synthetic and real data from Yahoo.

A cut-down version of this package under a GPL licence is available from http://github.com/robjhyndman/anomalous.

Simple Example
==============

````
  z <- ts(matrix(rnorm(3000),ncol=100),freq=4)
  y <- tsmeasures(z)
  biplot.features(y)
  anomaly(y)
````

ACM Licence
