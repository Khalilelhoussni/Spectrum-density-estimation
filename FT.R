FT=function(x) {
  n <- length(x)
  ct0=sqrt(2/n)
  seq0= c(1,n/2 + 1)
  ft <- fft(x)
  ft1 <- ft[1:(n/2 + 1)]
  
  c_values <- Im(ft1[-seq0])
  d_values <- Re(ft1[-seq0])
  
  Z <- numeric(n)
  Z[1] <- ft1[1]/sqrt(n)
  
  Z[seq(2, n-2, by = 2)] <- -ct0 * c_values
  Z[seq(3, n-1, by = 2)] <-  ct0 * d_values
  
  Z[n] <- ft1[n/2+1]/sqrt(n)
  
  return(Z)
}
