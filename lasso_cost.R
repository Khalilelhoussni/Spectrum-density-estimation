lasso_cost = function(X, y, beta, lambda, type="LASSO", intercept=T) {
  p=ncol(X)
  n=length(y)
  lbeta=length(beta)
  if((lbeta-intercept)==p){
    i0=as.numeric(intercept)
    #residuals <- X %*% beta[2:(2*n+1)] - y
    beta1= beta[i0+1:(2*n)]
    residuals <- X %*%beta1 - y
    if(intercept) residuals=residuals+beta[1]
    #cost <- 0.5*sum(residuals^2) + lambda * sum(abs(beta[2:(2*n+1)]))
    if(type=="LASSO") cost = 0.5*sum(residuals^2) + lambda * sum(abs(beta1))
    if(type=="sqrt-LASSO") cost = sqrt(sum(residuals^2)) + lambda * sum(abs(beta1))
  }
  else warning("Dimension problems in X, beta and intercept")
  return(cost)
}
