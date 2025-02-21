spectrum_plot <- function (z1,z2,z3,y,alpha,beta,freq0,coefAR,coefMA,sd=sd,type="Atomic",p0=p0,intercept=T, filter.number=filter.number, family=family, bc="periodic",L){
  i0 = as.numeric(intercept)
  n <- (length(z1)-i0)/2
  p <- 2*n
  yt0 <- wd(rep(0,n), filter.number=filter.number, family=family, bc="periodic")
  
  spectr4 <- numeric(n/2)
  spectr5 <- numeric(n/2)
  spectr6 <- numeric(n/2)
  spectr7 <- numeric(n/2)
  
  freq <- seq(1/n , 0.5 , length.out = n/2 )
  spe <- numeric(n/2)
  j <- 1
  if (type == "All"){
    beta_phi2 = wr_vec(z2[i0+((n+1):p)], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")
    beta_phi1 = wr_vec(z1[i0+((n+1):p)], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")
    beta_phi3 = wr_vec(z3[i0+((n+1):p)], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")
    
    spectr1 = (z1[i0+(1:n)] + beta_phi1)^2
    spectr2 = (z2[i0+(1:n)] + beta_phi2)^2
    spectr3 = (z3[i0+(1:n)] + beta_phi3)^2
    
    for (i in (1:(n/2)) ){
      spe[i] <- Spectrum_real(coefAR,coefMA,sd,freq[i])
      print(spe[i])
      
      if (freq[i] %in% freq0){
        indice = which(freq0 == freq[i])
        spe[i] = spe[i] + n*0.5*(alpha[indice]^2 + beta[indice]^2)
      }
    }
  }
  if (type == "Atomic"){
    spectr1 = (z1[i0+(1:n)])^2
    spectr2 = (z2[i0+(1:n)])^2
    spectr3 = (z3[i0+(1:n)])^2
    for (i in freq0){
      indice = which(freq0 == i)
      indice2 = which(freq == i)
      spe[indice2] = n*0.5*(alpha[indice]^2 + beta[indice]^2)
    }
  }
  if (type == "Smooth"){
    beta_phi2 = wr_vec(z2[i0+((n+1):p)], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")
    beta_phi1 = wr_vec(z1[i0+((n+1):p)], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")
    beta_phi3 = wr_vec(z3[i0+((n+1):p)], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")
    
    spectr1 = (beta_phi1)^2
    spectr2 = (beta_phi2)^2
    spectr3 = (beta_phi3)^2
    
    for (i in freq ){
      spe[j] <- Spectrum_real(coefAR,coefMA,sd,i)
      j<- j+1
    }
  }
  
  for (i in 1:((n/2)-1)) {
    spectr4[i] <- spectr1[2*i] + spectr1[2*i+1]
  }
 
  
  
  for (i in 1:((n/2)-1)) {
    spectr5[i] <- spectr2[2*i] + spectr2[2*i+1]
  }
  
  
  for (i in 1:((n/2)-1)) {
    spectr6[i] <- spectr3[2*i] + spectr3[2*i+1]
  }
  
  for (i in 1:L) {
    spectr7[i] <- spectr4[i]
  }
  
  
  for (i in (L+1):((n/2)-1)) {
    spectr7[i] <- spectr4[i]
    for (k in 1:L) {
      spectr7[i] <- spectr4[i-k] + spectr4[i+k]
    }
    spectr7[i] <- spectr7[i]/(2*L)
  }
  
  
  
  logYN='yes'
  spectre <- spectrum(y , log = logYN )
  lines(freq, spectr4, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "red")
  lines(freq,spe,type = "l",col= "orange")

         
  
  spectre <- spectrum(y , log = logYN )
  lines(freq, spectr6, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "green")
  lines(freq,spe,type = "l",col= "orange")

  
  spectre <- spectrum(y , log = logYN )
  lines(freq, spectr5, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "blue")
  lines(freq,spe,type = "l",col= "orange")
  
  spectre <- spectrum(y , log = logYN )
  lines(freq, spectr7, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "red")
  lines(freq,spe,type = "l",col= "orange")
  
  
}
