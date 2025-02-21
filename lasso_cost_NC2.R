lasso_cost_NC2 = function(X, y, beta, lambda, intercept=F,a =0.5) {
  p=ncol(X)
  n=length(y)
  lbeta=length(beta)
  i0 = as.numeric(intercept)
  if((lbeta-intercept)==p){
    i0=as.numeric(intercept)
    #residuals <- X %*% beta[2:(2*n+1)] - y
    beta1= beta[i0+(1:p)]
    residuals <- X %*%beta1 - y
    if(intercept) residuals=residuals+beta[1]
    #cost <- 0.5*sum(residuals^2) + lambda * sum(abs(beta[2:(2*n+1)]))
    cost = 0.5*sum(residuals^2) + lambda * sum(log(1 + abs(beta1)/a))
  }
  else warning("Dimension problems in X, beta and intercept")
  return(cost)
}
