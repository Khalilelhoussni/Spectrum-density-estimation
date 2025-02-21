lambda_QUT <- function(X,m,type,sd=NA,intercept=T,alpha=0.05) {
  ## Cacluate QUT for the cost function given in lasso_cost, either LASSO or sqrt-LASSO, but WITH an intercept
  ## X=input matrix

  out=NULL
  if (type == "LASSO" & is.na(sd)){
    warning("For LASSO, 'sd' must be specified")
  } 
  else {
    n=nrow(X)
    if (type == "LASSO"){ 
      quantiles_estimates <- numeric(m)
      z = matrix(rnorm(n*m,sd=sd), ncol=m)
      
      if (intercept){
        zm <- apply(z,2,mean)
        z  = t(t(z)-zm)
      }
      
      ty = t(X)%*%z
      nr = apply(abs(ty),2,max)
      quantiles_estimates = nr
    }
    
    if (type == "sqrt-LASSO") {
      quantiles_estimates <- numeric(m)
    
      z2 = matrix(rnorm(n*m,sd=1), ncol=m)
      
      if (intercept){
        zm2 <- apply(z2,2,mean)
        z2  = t(t(z2)-zm2)
      }
      ty2 = t(X)%*%z2
      nr2 = apply(abs(ty2),2,max)
      den=sqrt(apply(z2^2,2,sum))
      quantiles_estimates = t(t(nr2)/den)
    }
    
    lambdaQUT <- quantile(quantiles_estimates, 1-alpha)
  
    out$lambdas=quantiles_estimates
    out$lambdaQUT=lambdaQUT
    out$alpha=alpha
  }
  
  return(out)
  
}
