lasso_cost_fast = function( y, beta, lambda, type="LASSO", intercept=T, p0=4, filter.number=1) {
  n=length(y)
  p= 2*n
  lbeta=length(beta)
  if((lbeta-intercept)==p){
    i0=as.numeric(intercept)
    #residuals <- X %*% beta[2:(2*n+1)] - y
    yt1=wd(rep(0,n),filter.number=filter.number, family="DaubExPhase", bc="periodic")
    beta1= beta[i0+(1:p)]
    Zn <- FTinv(beta1[1:n])
    Zm <- Re(Zn)
    Rc <- FTinv(wr_vec(beta1[-c(1:n)], p0=p0, yt1, filter.number=filter.number, family="DaubExPhase", bc="periodic"))
    R <- Re(Rc)
    Xbeta1 <- Zm + R
    residuals <- Xbeta1 - y
    if(intercept) residuals=residuals+beta[1]
    #cost <- 0.5*sum(residuals^2) + lambda * sum(abs(beta[2:(2*n+1)]))
    if(type=="LASSO") cost = 0.5*sum(residuals^2) + lambda * sum(abs(beta1))
    if(type=="sqrt-LASSO") cost = sqrt(sum(residuals^2)) + lambda * sum(abs(beta1))
  }
  else warning("Dimension problems in X, beta and intercept")
  return(cost)
}
