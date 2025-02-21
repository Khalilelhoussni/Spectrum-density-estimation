Q_fista <- function(y,x,z,L,lambda,intercept=T, type = "LASSO") {
  p0=4
  filter.number=1
  J=log2(p0)
  n=length(y)
  p= 2*n
  lbeta=length(z)
  if((lbeta-intercept)==p){
    i0=as.numeric(intercept)
    #residuals <- X %*% beta[2:(2*n+1)] - y
    beta1 = z[i0+(1:p)]
    beta2 = x[i0+(1:p)]
    yt1 = wd(rep(0,n),filter.number=filter.number, family="DaubExPhase", bc="periodic")
    Zn <- FTinv(beta1[1:n])
    Zm <- Re(Zn)
    Rc <- FTinv(wr_vec(beta1[-c(1:n)], p0=p0, yt1, filter.number=filter.number, family="DaubExPhase", bc="periodic"))
    R <- Re(Rc)
    Xbeta1 <- Zm + R
    
    residuals <- Xbeta1 - y
    if(intercept) residuals=residuals+z[1]
    
    Zc2 <- Re(FT(residuals))
    Rc2 <- Re(wd_vec(Zc2, p0, filter.number, family="DaubExPhase", bc="periodic"))
    gradient1 <- c(Zc2,Rc2) / sqrt(sum((residuals)^2))
    gradient2 <- c(Zc2,Rc2) 
    
    if (type == "sqrt-LASSO") Q <- sqrt(sum(residuals^2)) + sum((beta2-beta1)*gradient1) + L/2 * sum((beta2-beta1)*(beta2-beta1)) + lambda * sum(abs(beta2))
    if (type == "LASSO")  Q <- 0.5*sum(residuals^2) + sum((beta2-beta1)*gradient2) + L/2 * sum((beta2-beta1)*(beta2-beta1)) + lambda * sum(abs(beta2))
  }
  else warning("Dimension problems in X, beta and intercept")
  return(Q)
  
}
 