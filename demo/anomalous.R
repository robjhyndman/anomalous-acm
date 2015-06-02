# Compute feature matrix for the real "dat0" data set
features0 <- tsmeasures(dat0, width = 24, window = 48)
# biplot on robust PCA for the "feature" matrix
biplot(features0)
# Anomaly detection by bivariate kernel density
hdr <- anomaly(features0, n = 3, plot = TRUE, method = "hdr")
# Anomaly detection by bivariate alpha-hull
ahull <- anomaly(features0, n = 3, plot = TRUE, method = "ahull")
# Plot top 3 anomalous time series by bivariate kernel density
plot(dat0[, hdr$index])
# Plot top 3 anomalous time series by ahull 
plot(dat0[, ahull$index])

# Compute feature matrix for the synthetic "dat4" data set
features0 <- tsmeasures(dat4, width = 24, window = 48)
# biplot on robust PCA for the "feature" matrix
biplot(features0)
# Anomaly detection by bivariate kernel density
hdr <- anomaly(features0, n = 3, plot = TRUE, method = "hdr")
# Anomaly detection by bivariate alpha-hull
ahull <- anomaly(features0, n = 3, plot = TRUE, method = "ahull")
# Plot top 3 anomalous time series by bivariate kernel density
plot(dat4[, hdr$index])
# Plot top 3 anomalous time series by ahull 
plot(dat4[, ahull$index])

# Compute feature matrix for the synthetic "dat5" data set
features0 <- tsmeasures(dat5, width = 24, window = 48)
# biplot on robust PCA for the "feature" matrix
biplot(features0)
# Anomaly detection by bivariate kernel density
hdr <- anomaly(features0, n = 3, plot = TRUE, method = "hdr")
# Anomaly detection by bivariate alpha-hull
ahull <- anomaly(features0, n = 3, plot = TRUE, method = "ahull")
# Plot top 3 anomalous time series by bivariate kernel density
plot(dat5[, hdr$index])
# Plot top 3 anomalous time series by ahull 
plot(dat5[, ahull$index])
