projsplx <- function(y) {
  n <- length(y)
  s <- sort(y, decreasing = TRUE)
  tmpsum <- cumsum(s[1:(n - 1)])
  tmp <- (tmpsum - 1) / (1:(n - 1))
  ind <- which(tmp >= s[2:n])[1]
  if (!is.na(ind)) {
    t <- tmp[ind]
  } else {
    t <- (tmpsum[n - 1] + s[n] - 1) / n
  }
  x <- apply(cbind(y - t, 0), 1, max)
  return(x)
}
