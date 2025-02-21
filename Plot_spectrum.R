
yt0=wd(rep(0,n),filter.number=filter.number, family="DaubExPhase", bc="periodic")

# computation of 
beta_phi2 = wr.vec(testFISTA1$beta[((n+1):length(outfit4$beta))], p0=p0, yt0, filter.number=filter.number, family="DaubExPhase", bc="periodic")
beta_phi1 = wr.vec(outfit3$beta[((n+1):length(outfit3$beta))], p0=p0, yt0, filter.number=filter.number, family="DaubExPhase", bc="periodic")


spectr = 0.5*(outfit3$beta[1:n] + beta_phi1)^2
spectr2 = 0.5*(outfit4$beta[1:n] + beta_phi2)^2
# Initialisez un vecteur pour stocker la somme des paires d'éléments consécutifs
spectr1 <- numeric(length(spectr) %/% 2)
spectr3 <- numeric(length(spectr2) %/% 2)


# Sommez chaque paire d'éléments consécutifs
#spectr1[1] <- 2*spectr[1]
for (i in 1:((length(spectr) %/% 2)-1)) {
  spectr1[i] <- spectr[2*i] + spectr[2*i+1]
}
spectr1[length(spectr) %/% 2] <- 2*spectr[length(spectr) %/% 2]


for (i in 1:((length(spectr2) %/% 2)-1)) {
  spectr3[i] <- spectr2[2*i] + spectr2[2*i+1]
}
spectr3[length(spectr2) %/% 2] <- 2*spectr[length(spectr2) %/% 2]


spectr3[1] <- 2*spectr2[1]


# Fréquences associées au spectre
freq <- seq(1/n , 0.5 , length.out = n/2 )

# Tracé du spectre

spectre <- spectrum(y , log = 'no' )


#spectree <-  (0.5*n*(0.05^2)) / abs(1 - 0.9 * exp(-2 * pi * 1i * freq))^2 + atomic

#plot(freq, spectree, type = "l", xlab = "Fréquence", ylab = "Densité Spectrale de Puissance", main = "Densité Spectrale ")


lines(freq, spectr1, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "red" )
lines(freq, spectr3, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "blue" )


