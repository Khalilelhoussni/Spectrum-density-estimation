FTinv <- function(x) {
  n <- length(x)
  yc1 <- numeric(n/2)
  yc2 <-  numeric(n/2)
  yc3 <- numeric(n/2)
  yc4 <-  numeric(n/2)
  yc1[2:(n/2)] <- x[seq(2, n-2, by = 2)]
  yc2 <- c(x[1], x[seq(3, n-1, by = 2)])
  yc3 <- numeric(n/2)
  yc4 <- numeric(n/2)
  
  ys <- c(yc1, yc3)
  yc <- c(yc2, yc4)
  
  ft1 <- fft(ys)
  ft2 <- fft(yc)
  
  k <- 1:n
  
  c_values <- -Im(ft1[k])
  d_values <- Re(ft2[k])
  
  Z <- ((1-sqrt(2))/sqrt(n)) * x[1] + sqrt(2/n) * (d_values + c_values) + (cos(pi*(k-1))/sqrt(n)) * x[n]
  
  return(Z)
}
