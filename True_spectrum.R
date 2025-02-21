Spectrum_real <- function(coefAR,coefMA,sd,frequence) {
  
  omega <-  -1i*2*pi*frequence
  p <- length(coefAR)
  q <- length(coefMA)
  polyAR <- 1
  polyMA <- 1
  
  for (i in 1:p) {
    polyAR <- polyAR - coefAR[i]*exp(omega)^i
  }
  for (i in 1:q) {
    polyMA <- polyMA + coefMA[i]*exp(omega)^i
  }
  #polyAR <- 1 - sum(sapply(1:p, function(i) coefAR[i] * exp(omega)^i))
  #polyMA <- 1 + sum(sapply(1:q, function(i) coefMA[i] * exp(omega)^i))
  
  density <- (sd^2/(pi)) * (abs(polyMA)/abs(polyAR))^2
  return(density)
  
}


     