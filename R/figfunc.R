biplot.features <- function(x, robust = TRUE, col, ...) {
  nc <- nrow(x)
  naomit.x <- na.omit(x)
  na.act <- na.action(naomit.x)
  if (is.null(na.act)) {
    avl <- 1:nc
  } else {
    avl <- (1:nc)[-na.action(naomit.x)]
  }
  if (missing(col)) {
    col <- c("#000000", "darkred")
  } else {
    lencol <- length(col)
    if (lencol == 1L) {
      col <- rep(col, 2)
    } else {
      col <- unique(col)[1:2]
    }
  }
  if (robust) {
    rbt.pca <- pcaPP::PCAproj(naomit.x, k = 2, scale = sd, center = mean)
  } else {
    rbt.pca <- princomp(naomit.x, cor = TRUE)
  }
  biplot(rbt.pca, col = col, xlabs = avl, xlab = "PC1", ylab = "PC2", 
         cex = c(0.8, 1), ...)
}
