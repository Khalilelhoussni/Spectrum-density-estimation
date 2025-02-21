rm(list=ls())
setwd("C:/Users/Khalil/Dropbox/ProjetsStatMachLearning 2023-24/Spectrum")
setwd("C:/Users/sylvain/Dropbox/Teaching/ProjetsStatMachLearning 2023-24/Spectrum")


library(glmnet)
library("wavethresh")


source("lambda_QUT.R")
source("lambda_QUT_fast.R")

source("lasso_cost.R")
source("lasso_cost_fast.R")

source("soft_thresholding.R")
source("FISTA_LASSO_SQRT-LASSO.R")
#source("FISTA_LASSO_SQRT-LASSO_fast.R")
source("FISTA_LASSO_SQRT-LASSO_fast.R")
source("ISTA.R")
source("FT.R")
source("FTinv.R")
source("wd_vec.R")
source("wr_vec.R")
source("BCR-lasso.R")
source("BCRfast-lasso.R")
source("Q_fista.R")

source("fista_nonconvex.R")
source("fista_nonconvex2.R")
source("FISTA_Linesearch2.R")
source("Thresholding_function_NC.R")
source("Thresholding_function_NONconvex.R")
source("lasso_cost_NC.R")
source("lasso_cost_NC2.R")

source("True_spectrum.R")
source("Plot_spectre.R")




n=1024
freq0s=c(0.1015625000,0.2656250000,0.3984375000)
alphas=c(3,2,3)/2
betas=c(2,3,2)/2

x <- c(1:n)
sd=1

mu=0
mu=arima.sim(n = n, list(ar = c(-0.2,-0.5,-0.1),ma = c(0.2,0.5,0.2,-0.5), sd = sd))

mu=as.vector(mu)

for(j in 1:3){
  mu=mu+alphas[j]*cos(2*pi*freq0s[j]*(x-1))+betas[j]*sin(2*pi*freq0s[j]*(x-1))
}

y=mu 

Xc <- matrix(0, nrow = n, ncol = n)
Xc[, 1] <- 1 / sqrt(n)
for (k in 1:(n/2 - 1)) {
  Xc[, 2*k] <- sqrt(2/n) * sin(2 * pi * (x-1) * (k / n))
  Xc[, 2*k + 1] <- sqrt(2/n) * cos(2 * pi * (x-1) * (k / n))
}
Xc[, n] <- sqrt(1/n) * cos(pi*(x-1))


p0=4
filter.number=1
family = "DaubExPhase"

Xw=matrix(NA, n, n)
J=log2(p0)
for(i in 1:n){
  yX=wd_vec(Xc[i,], p0, filter.number=1, family=family, bc="periodic")
  Xw[,i]=yX
}
Xphi = t(Xw)

X = cbind(Xc,Xphi)

m=1000
outQUT=lambda_QUT(X,m,type="LASSO",intercept=T,sd=sd,alpha=0.05) 
lambdaQUT=outQUT$lambdaQUT
outQUTlasso=lambda_QUT_fast(n,m,type="LASSO", sd = sd,intercept=T, alpha=0.05,p0=p0,filter.number=1,family=family) 
lambdaQUTlasso=outQUTlasso$lambdaQUT

outQUTsqrtlasso=lambda_QUT_fast(n,m,type="sqrt-LASSO",intercept=F,sd=sd,alpha=0.05) 
lambdaQUTsqrtlasso=outQUTsqrtlasso$lambdaQUT

tol=1.e-16

beta=rnorm(2*n+1)
lasso_cost_fast( y, beta, lambdaQUTlasso, type="LASSO", intercept=T, p0=p0, filter.number=filter.number) 
lasso_cost_NC(X, y, beta, lambdaQUTlasso, intercept=T,eps=1)

testBCR1=BCRfast_lasso(y, lambda=lambdaQUTlasso, intercept=T, savecost=T, tol = 1e-20)

testFISTA1=fista_lasso_sqrt_lasso_fast(y, lambda=lambdaQUTlasso/6, intercept=T, type ="LASSO", tol = tol, rho=1.2, L0=0.1)
testFISTA2=fista_lasso_sqrt_lasso_fast4(y, lambda=lambdaQUTlasso, intercept=T, type ="LASSO", tol = tol, rho=1.1, L0=1)
testFISTA3=fista_lasso_sqrt_lasso_fast4(y, lambda=lambdaQUT2, intercept=T, type ="sqrt-LASSO", tol = tol, rho=1.1, L0=0.1)
testFISTA4=fista_lasso_sqrt_lasso_fast3(y, lambda=lambdaQUT2, intercept=T, type ="sqrt-LASSO", tol = tol, rho=1.1, L0=1000)
lambda0=max(abs(t(X)%*%(y-mean(y))))
plot(testFISTA1$beta)
plot(Sol8bis$beta-Sol10$beta)
plot(testFISTA1$beta-testFISTA2$beta)
plot(Sol5$beta-testFISTA3$beta)

# LASSO WITH MATRIX
Sol1 <- fista_lasso2(X,y,lambda=lambdaQUT, intercept =F, tol=tol)

Sol1bis=BCR_lasso_test(Xc, Xphi, y, lambda=lambdaQUT, intercept=T, tol = 1e-20)

SOL <- fista_lasso_sqrt_lasso(X, y, m=1000,sd=1, intercept=T,type ="sqrt-LASSO" ,tol = tol)
SOL2 <- fista_lasso_sqrt_lasso_fast(y, intercept=T,type ="LASSO" ,tol = tol)

# LASSO WITH FFT AND DWT
Sol2 <- fista_lasso4(y,lambda=lambdaQUT, tol=tol)
# LASSO  BCR WITH MATRIX
Sol3 <-  BCR_lasso(Xc,Xphi, y, lambda=lambdaQUT, tol=tol)
plot(log(Sol3$cost), type="l")
# LASSO BCR WITH FFT AND DWT
Sol4 <- BCR_lasso_fft(y, lambda=lambdaQUT, tol=tol)
# sqrt_LASSO WITH MATRIX
Sol5 <-  fista_sqlasso(X, y, lambda = lambdaQUT2  ,intercept =T, tol=1.e-20)
# sqrt_LASSO WITH FFT AND DWT
Sol6 <- fista_sqlasso_fft(y, lambda=lambdaQUT2, tol=tol)
# LASSO WITH GLMNET
Sol7 <- glmnet(X, y, alpha = 1, standardize = F, intercept = T,  lambda = 2*(4.05)/length(y), thresh=1.e-20)
# LASSO/LASSO_BCR WITH intercept
Sol8 <- fista_lasso_intercept(X, y, lambda= lambdaQUT, tol=tol)

Sol8bis <- fista_lasso_intercept_test(X, y, lambdaQUT , intercept=T, tol = tol)

Sol9 <- BCR_lasso_intercept(Xc,Xphi, y, lambda=lambdaQUT, tol=tol)

Sol10 <- fista_nonconvex(X, y, lambdaQUTlasso/6, intercept=T, tol = tol,a=2)
Sol11 <- fista_nonconvex2(X, y, lambdaQUTlasso/6, intercept=T, tol = tol,nu=0.1)


plot(Sol10$beta)
lines(Sol11$beta, col = "red")
lines(testFISTA1$beta, col = "blue")
# Comparison between the solutions of  LASSO_FISTA and LASSO_FISTA with FFT and DWT

plot(Sol1$beta-Sol2$beta)

# Comparison between the solutions LASSO_BCR and LASSO_BCR with FFT and DWT

plot(Sol2$beta-Sol4$beta)

# Comparison between the solutions of SQRT_LASSO_FISTA and SQRT_LASSO_FISTA with FFT and DWT 

plot(Sol5$beta-SOL$beta)
plot(Sol5$beta)

# Comparison between the solutions of SQRT_LASSO_FISTA and LASSO_FISTA ( value of L and lambdaqut to verify ! )

plot(Sol5$beta-Sol1$beta)

# Comparison between the solutions of LASSO_FISTA and LASSO_BCR 

plot(Sol1$beta-Sol3$beta)

# Comparison between the solutions of LASSO_FISTA and LASSO_BCR with intercept

plot(Sol8bis$beta-SOL$beta)

# Comparison between the solutions of glmnet and LASSO_FISTA 

plot(Sol8bis$beta[2:(2*n+1)]-Sol7$beta)
plot(Sol9$beta[2:(2*n+1)])




# Spectrum estimation plot

testFISTA2=fista_lasso_sqrt_lasso_fast(y, lambda=lambdaQUTlasso/15, intercept=F, type ="LASSO", tol = tol, rho=2.0, L0=0.1)
testFISTA3=fista_lasso_sqrt_lasso_fast(y, lambda=lambdaQUTsqrtlasso/15, intercept=F, type ="sqrt-LASSO", tol = tol, rho=2.0, L0=0.1)

yt0 <- wd(rep(0,n), filter.number=filter.number, family=family, bc="periodic")
beta_phi2 = wr_vec(testFISTA1$beta[((n+1):length(testFISTA1$beta))], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")
beta_phi1 = wr_vec(Sol10$beta[((n+1):length(Sol10$beta))], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")
beta_phi3 = wr_vec(Sol11$beta[((n+1):length(Sol11$beta))], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")


spectr1 = (4/n)*(Sol10$beta[1:n] + beta_phi1)^2
spectr2 = (4/n)*(testFISTA1$beta[1:n] + beta_phi2)^2
spectr3 = (4/n)*(Sol11$beta[1:n] + beta_phi3)^2

spectr4 <- numeric(length(spectr1) %/% 2)
spectr5 <- numeric(length(spectr2) %/% 2)
spectr6 <- numeric(length(spectr3) %/% 2)


for (i in 1:((length(spectr1) %/% 2)-1)) {
  spectr4[i] <- spectr1[2*i] + spectr1[2*i+1]
}
spectr4[length(spectr1) %/% 2] <- 2*spectr1[length(spectr1) %/% 2]


for (i in 1:((length(spectr2) %/% 2)-1)) {
  spectr5[i] <- spectr2[2*i] + spectr2[2*i+1]
}
spectr5[length(spectr2) %/% 2] <- 2*spectr2[length(spectr2) %/% 2]

for (i in 1:((length(spectr3) %/% 2)-1)) {
  spectr6[i] <- spectr3[2*i] + spectr3[2*i+1]
}
spectr6[length(spectr3) %/% 2] <- 2*spectr3[length(spectr3) %/% 2]




# Fréquences associées au spectre
freq <- seq(1/n , 0.5 , length.out = n/2 )

# Tracé du spectre

par(mfrow=c(1,4))
logYN='no'
#spectre <- spectrum(y , log = logYN )
#spectree <-  (0.5*n*(0.05^2)) / abs(1 - 0.9 * exp(-2 * pi * 1i * freq))^2 + atomic
#plot(freq, spectree, type = "l", xlab = "Fréquence", ylab = "Densité Spectrale de Puissance", main = "Densité Spectrale ")
spe <- numeric(512)
j <- 1
for (i in freq ){
  spe[j] <- Spectrum_real( c(-0.2,-0.5,-0.1), c(0.2,0.5,0.2,-0.5),1,i)
  j<- j+1
}
plot(freq,spe)
plot(freq, spectr4, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "red")
plot(freq, spectr6, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "green" )
plot(freq, spectr5, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "blue" )


par(mfrow=c(1,1))
#spectre <- spectrum(y , log = logYN )

spe <- numeric(512)
j <- 1
for (i in freq ){
  spe[j] <- Spectrum_real( c(-0.2,-0.5,-0.1), c(0.2,0.5,0.2,-0.5),1,i)
  j<- j+1
}
plot(freq,spe)

lines(freq, spectr4, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "red")
lines(freq, spectr6, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "green")
lines(freq, spectr5, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "blue")

xt <- testFISTA1$beta
yt <- Sol10$beta
zt <- Sol11$beta

a3 <- c(-0.2,-0.5,-0.1)
a4 <- c(0.2,0.5,0.2,-0.5)

spectrum_plot(xt,yt,zt,y,alphas,betas,freq0s,a3,a4,sd=sd,type="All",p0=p0,intercept=T, filter.number=filter.number, family=family, bc="periodic",L=10)

# Non_equispaced time serie 

proportion_NA <- 0.1
yn <- y
indices_missing <- sort(sample(1:n, round(n*proportion_NA)))
indices_missing <- as.numeric(indices_missing)
indices_value <- setdiff(1:n,indices_missing)
indices_value <- as.numeric(indices_value)

yn[indices_missing] <- NA

yn_cleaned <- yn[!is.na(yn)]
Ln <- length(yn_cleaned)
nNA <- 2^floor(log(Ln,2))
yn2 <- yn_cleaned[1:nNA]
indices_value2 <- indices_value[1:nNA]

X_NA <- X[indices_value,]


Sol1_NA <- fista_lasso2(X_NA,yn_cleaned,intercept = F, lambda=lambdaQUT/2, tol=tol)
Sol1 <- fista_lasso2(X,y,lambda=lambdaQUT/2, intercept =F, tol=tol)


SS <- X%*%Sol1_NA$beta
yn[indices_missing] <- SS[indices_missing]
print(y[indices_missing]-yn[indices_missing])
print(y[indices_missing])

# Spectrum missing data vs full data

yt0=wd(rep(0,n),filter.number=filter.number, family=family, bc="periodic")



beta_phi2bis = wr_vec(Sol1$beta[((n+1):length(Sol1$beta))], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")
beta_phi1bis = wr_vec(Sol1_NA$beta[((n+1):length(Sol1_NA$beta))], p0=p0, yt0, filter.number=filter.number, family=family, bc="periodic")


spectrbis = 0.5*(Sol1_NA$beta[1:n] + beta_phi1bis)^2
spectr2bis = 0.5*(Sol1$beta[1:n] + beta_phi2bis)^2
# Initialisez un vecteur pour stocker la somme des paires d'éléments consécutifs
spectr1bis <- numeric(length(spectrbis) %/% 2)
spectr3bis <- numeric(length(spectr2bis) %/% 2)


# Sommez chaque paire d'éléments consécutifs
#spectr1[1] <- 2*spectr[1]
for (i in 1:((length(spectrbis) %/% 2)-1)) {
  spectr1bis[i] <- spectrbis[2*i] + spectrbis[2*i+1]
}
spectr1bis[length(spectrbis) %/% 2] <- 2*spectrbis[length(spectrbis) %/% 2]


for (i in 1:((length(spectr2bis) %/% 2)-1)) {
  spectr3bis[i] <- spectr2bis[2*i] + spectr2bis[2*i+1]
}
spectr3bis[length(spectr2bis) %/% 2] <- 2*spectrbis[length(spectr2bis) %/% 2]


# Fréquences associées au spectre
freq <- seq(1/n , 0.5 , length.out = n/2 )

# Tracé du spectre

spectrebis <- spectrum(y , log = 'no' )


#spectree <-  (0.5*n*(0.05^2)) / abs(1 - 0.9 * exp(-2 * pi * 1i * freq))^2 + atomic

#plot(freq, spectree, type = "l", xlab = "Fréquence", ylab = "Densité Spectrale de Puissance", main = "Densité Spectrale ")


lines(freq, spectr1bis, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "red" )
lines(freq, spectr3bis, type = "l", xlab = "Fréquence", ylab = "Spectre", col = "blue" )


