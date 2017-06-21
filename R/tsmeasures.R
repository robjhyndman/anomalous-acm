tsmeasures <- function(y, normalise = TRUE, 
  width=ifelse(frequency(y)>1, frequency(y), 10), window=width) {
  # y: a multivariate time series
  # normalise: TRUE: scale data to be normally distributed
  # width: a window size for variance change and level shift, lumpiness
  # window: a window size for KLscore
  y <- as.ts(y)
  tspy <- tsp(y)
  freq <- frequency(y)
  if (width <= 1L | window <= 1L)
    stop("window widths should be greater than 1.")

  #if single time series is supplied, convert to matrix
  if(is.null(dim(y))) {
	y <- as.ts(as.matrix(y))
	tsp(y) <- tspy
  } 

  # Remove columns containing all NAs
  nay <- is.na(y)
  allna <- apply(nay, 2, all)
  x <- y[, !allna, drop=FALSE]
  if (normalise) {
    # Normalise data to mean = 0, sd = 1, unless constant, then normalize to all 0
    x <- as.ts(apply(x, 2, function(x) {
      scale(x, center = TRUE, scale = !is.constant(x))
    }))

  }
  trimx <- as.ts(apply(x, 2, Trim))
  tsp(trimx) <- tsp(x) <- tspy
  measures <- list()
  measures$lumpiness <- apply(x, 2, Lump, width = width)
  measures$entropy <- apply(x, 2, Entropy)
  measures$ACF1 <- apply(x, 2, FoAcf)
  measures$lshift <- apply(trimx, 2, RLshift, width = width)
  measures$vchange <- apply(trimx, 2, RVarChange, width = width)
  measures$cpoints <- apply(x, 2, Cpoints)
  measures$fspots <- apply(x, 2, Fspots)
#  measures$mean <- colMeans(x, na.rm = TRUE)
#  measures$var <- apply(x, 2, var, na.rm = TRUE)
  varts <- apply(x, 2, VarTS, tspx = tspy)
  measures$trend <- sapply(varts, function(x) x$trend)
  measures$linearity <- sapply(varts, function(x) x$linearity)
  measures$curvature <- sapply(varts, function(x) x$curvature)
  measures$spikiness <- sapply(varts, function(x) x$spike)
  if (freq > 1) {
    measures$season <- sapply(varts, function(x) x$season)
    measures$peak <- sapply(varts, function(x) x$peak)
    measures$trough <- sapply(varts, function(x) x$trough)
  }
  if (nrow(y) <= (2 * window)) {
    #warning("I cannot compute KLscore when the length is too small.")
		measures$KLscore <- NA
		measures$change.idx <- NA
		
  } else {
    threshold <- dnorm(38)
    kl <- apply(x, 2, function(x) 
                KLscore(x, window = window, threshold = threshold))
    measures$KLscore <- sapply(kl, function(x) x$score)
    measures$change.idx <- sapply(kl, function(x) x$change.idx)
  }
  tmp <- do.call(cbind, measures)
  nr <- ncol(y)
  nc <- length(measures)
  mat <- matrix(, nrow = nr, ncol = nc)
  colnames(mat) <- colnames(tmp)
  mat[!allna, ] <- tmp
  out <- structure(mat, class = c("features", "matrix"))
  return(out)
}

#this function is copied from the forecast package
is.constant <- function (x) 
{
    x <- as.numeric(x)
    y <- rep(x[1], length(x))
    return(identical(x, y))
}

# Trimmed time series elimating outliers's influence
Trim <- function(x, trim = 0.1) {
  qtl <- quantile(x, c(trim, 1 - trim), na.rm = TRUE)
  lo <- qtl[1L]
  hi <- qtl[2L]
  x[x < lo | x > hi] <- NA
  return(x)
}

# Lumpiness: cannot be used for yearly data
Lump <- function(x, width) {
  start <- seq(1, nr <- length(x), by = width)
  end <- seq(width, nr + width, by = width)
  nsegs <- nr/width
  varx <- sapply(1:nsegs, function(idx)
                 var(x[start[idx]:end[idx]], na.rm = TRUE))
  lumpiness <- var(varx, na.rm = TRUE)
  return(lumpiness)
}

# spectral_entropy from ForeCA package
Entropy <- function(x) {
  entropy <- try(ForeCA::spectral_entropy(na.contiguous(x))[1L], silent = TRUE)
  if (class(entropy) == "try-error") {
    entropy <- NA
  }
  return(entropy)
}

# First order of autocorrelation
FoAcf <- function(x) {
  foacf <- try(acf(x, plot = FALSE, na.action = na.exclude)$acf[2L], silent = TRUE)
  if (class(foacf) == "try-error") {
    foacf <- NA
  }
  return(foacf)
}

# Level shift using rolling window
RLshift <- function(x, width) {
  rollmean <- try(RcppRoll::roll_mean(x, width, na.rm = TRUE), silent = TRUE)
  
  if (class(rollmean) == "try-error") {
    lshifts <- NA
  } else {
    lshifts <- tryCatch(max(abs(diff(rollmean, width)), na.rm = TRUE),
        warning = function(w) w)
    if (any(class(lshifts) == "warning")) {
      lshifts <- NA
    }
  }
  return(lshifts)
}

RVarChange <- function(x, width) {
  rollvar <- try(RcppRoll::roll_var(x, width, na.rm = TRUE), silent = TRUE)
  
  if (class(rollvar) == "try-error") {
    vchange <- NA
  } else {
    vchange <- tryCatch(max(abs(diff(rollvar, width)), na.rm = TRUE),
        warning = function(w) w)
    if (any(class(vchange) == "warning")) {
      vchange <- NA
    }
  }
  return(vchange)
}

# the number of crossing points
Cpoints <- function(x) {
  midline <- sum(range(x, na.rm = TRUE))/2
  ab <- x <= midline
  lenx <- length(x)
  p1 <- ab[1:(lenx - 1)]
  p2 <- ab[2:lenx]
  cross <- (p1 & !p2) | (p2 & !p1)
  return(sum(cross))
}

# Flat spots using disretization
Fspots <- function(x) {
  cutx <- try(cut(x, breaks = 10, include.lowest = TRUE, labels = FALSE),
              silent = TRUE)
  if (class(cutx) == "try-error") {
    fspots <- NA
  } else {
    rlex <- rle(cutx)
    # Any flat spot
    return(max(rlex$lengths))
    # Low flat spots
    # ones <- (rlex$values == 1)
    # return(max(rlex$lengths[ones]))
  }
}

# Strength of trend and seasonality and spike
VarTS <- function(x, tspx) {
  x <- as.ts(x)
  tsp(x) <- tspx
  freq <- tspx[3]
  contx <- try(na.contiguous(x), silent = TRUE)
  len.contx <- length(contx)
  if (length(contx) <= 2 * freq || class(contx) == "try-error") {
    trend <- linearity <- curvature <- season <- spike <- peak <- trough <- NA
  } else {
    if (freq > 1L) {
      all.stl <- stl(contx, s.window = "periodic", robust = TRUE)
      starty <- start(contx)[2L]
      pk <- (starty + which.max(all.stl$time.series[, "seasonal"]) - 1L) %% freq
      th <- (starty + which.min(all.stl$time.series[, "seasonal"]) - 1L) %% freq
      pk <- ifelse(pk == 0, freq, pk)
      th <- ifelse(th == 0, freq, th)
      trend0 <- all.stl$time.series[, "trend"]
      fits <- trend0 + all.stl$time.series[, "seasonal"]
      adj.x <- contx - fits
      v.adj <- var(adj.x, na.rm = TRUE)
      detrend <- contx - trend0
      deseason <- contx - all.stl$time.series[, "seasonal"]
      peak <- pk * max(all.stl$time.series[, "seasonal"], na.rm = TRUE)
      trough <- th * min(all.stl$time.series[, "seasonal"], na.rm = TRUE)
      remainder <- all.stl$time.series[, "remainder"]
      season <- ifelse(var(detrend, na.rm = TRUE) < 1e-10, 0,
                       max(0, min(1, 1 - v.adj/var(detrend, na.rm = TRUE))))
    } else { # No seasonal component
      tt <- 1:len.contx
      trend0 <- fitted(mgcv::gam(contx ~ s(tt)))
      remainder <- contx - trend0
      deseason <- contx - trend0
      v.adj <- var(remainder, na.rm = TRUE)
    }
    trend <- ifelse(var(deseason, na.rm = TRUE) < 1e-10, 0,
                    max(0, min(1, 1 - v.adj/var(deseason, na.rm = TRUE))))
    n <- length(remainder)
    v <- var(remainder, na.rm = TRUE)
    d <- (remainder - mean(remainder, na.rm = TRUE))^2
    varloo <- (v * (n - 1) - d)/(n - 2)
    spike <- var(varloo, na.rm = TRUE)
    pl <- poly(1:len.contx, degree = 2)
    tren.coef <- coef(lm(trend0 ~ pl))[2:3]
    linearity <- tren.coef[1]
    curvature <- tren.coef[2]
  }
  if (freq > 1) { 
    return(list(trend = trend, season = season, spike = spike,
                peak = peak, trough = trough, linearity = linearity,
                curvature = curvature))
  } else { # No seasonal component
    return(list(trend = trend, spike = spike, linearity = linearity,
                curvature = curvature))
  }
}

# Kullback-Leibler score
KLscore <- function(x, window, threshold = dnorm(38)) {
  gw <- 100 # grid width
  xgrid <- seq(min(x, na.rm = TRUE), max(x, na.rm = TRUE), length = gw)
  grid <- xgrid[2L] - xgrid[1L]
  tmpx <- x[!is.na(x)] # Remove NA to calculate bw
  bw <- bw.nrd0(tmpx) 
  lenx <- length(x)
  # Using binning algorithm to achieve efficiency but obsecure exact positions.
  # lastrep <- ceiling(lenx/5) 
  # group <- rep(1:lastrep, each = 5)[1:lenx]
  # midpoints <- aggregate(x, by = list(group), function(y) y[3L])[, 2]
  # dens.mat <- matrix(, nrow = lastrep, ncol = gw)
  # for (i in 1L:lastrep) {
  #   dens.mat[i, ] <- dnorm(xgrid, mean = midpoints[i], sd = bw)
  # }
  dens.mat <- matrix(, nrow = lenx, ncol = gw)
  for (i in 1L:lenx) {
    dens.mat[i, ] <- dnorm(xgrid, mean = x[i], sd = bw)
  }
  dens.mat <- pmax(dens.mat, threshold)
  rmean <- RcppRoll::roll_mean(dens.mat, n = window, na.rm = TRUE, fill = NA,
                               align = "right") # by column
  # lo <- seq(1, lastrep - window + 1)
  # hi <- seq(window + 1, lastrep)
  lo <- seq(1, lenx - window + 1)
  hi <- seq(window + 1, lenx)
  seqidx <- min(length(lo), length(hi))
  kl <- sapply(1:seqidx, function(i) sum(rmean[lo[i], ] *
               (log(rmean[lo[i], ]) - log(rmean[hi[i], ])) *
               grid, na.rm = TRUE))
  diffkl <- diff(kl, na.rm = TRUE)
  maxidx <- which.max(diffkl) + 1
  return(list(score = max(diffkl), change.idx = maxidx))
}
