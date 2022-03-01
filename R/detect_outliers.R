anomaly <- function(x, n = 10, method = c("hdr", "ahull"), robust = TRUE, 
    ordered=FALSE, plot = TRUE, labels = TRUE, col) {
  # x: a matrix returned by `tsmeasures` function
  nc <- nrow(x)
  if (nc < n) {
    stop("Your n is too large.")
  }
  x[is.infinite(x)] <- NA # ignore inf values
  naomit.x <- na.omit(x) # ignore missing values
  na.act <- na.action(naomit.x)
  if (is.null(na.act)) {
    avl <- 1:nc
  } else {
    avl <- (1:nc)[-na.action(naomit.x)]
  }
  method <- match.arg(method)
  # robust PCA space (scaling version)
  if (robust) {
    rbt.pca <- pcaPP::PCAproj(naomit.x, k = 2, center = mean, scale = sd)
  } else {
    rbt.pca <- princomp(scale(naomit.x, center = TRUE, scale = TRUE), 
                        cor = TRUE)
  }
  scores <- rbt.pca$scores[,1:2] #make non-robust PCA work, prevent dimension mismatch below
}

# Function to find the next nfind most outlying points in 2d pc space
# given those already found are indexed by outliers and
# the value of alpha corresponding to the last outlier was highalpha

findnextoutlier <- function(scores, highalpha, outliers, nfind=1)
{
  niter <- 500
  first <- 0
  last <- highalpha
  len.out <- 0
  numiter <- 0
  n <- length(outliers)+nfind
  while (len.out != n && (numiter <- numiter + 1) <= niter) {
    fit <- alphahull::ahull(scores, alpha = half <- (first + last)/2)
    radius <- fit$arcs[, 3]
    check <- radius == 0
    len.out <- length(radius[check])
    if (len.out >= 1 && len.out <= n) {
      xpos <- fit$arcs[check, 1]
      xidx <- which(is.element(scores[, 1], xpos))
      tmp.idx <- xidx
    }
    if (len.out >= 0 && len.out <= n) {
      last <- half
    } else {
      first <- half
    } 
  }
  if(numiter > niter)
    stop("Too hard. Please reduce the number of outliers requested")
  newoutlier <- setdiff(tmp.idx, outliers)
  if(length(newoutlier) != nfind)
    stop("Found too many or too few outliers")
  return(list(outlier=newoutlier, alpha=half))
}

