lambda_QUT_fast <- function(n,m,type,sd=NA,intercept = T,alpha=0.05,p0=1, filter.number=1, family="DaubExPhase", bc="periodic") {
  ## Cacluate QUT for the cost function given in lasso_cost, either LASSO or sqrt-LASSO, but WITH an intercept
  ## X=input matrix
  
  out=NULL
  if (type == "LASSO" & is.na(sd)){
    warning("For LASSO, 'sd' must be specified")
  } 
  else {
    if (type == "LASSO"){ 
      quantiles_estimates <- numeric(m)
      z = matrix(rnorm(n*m,sd=sd), ncol=m)
      if (intercept){
        zm <- apply(z,2,mean)
        z  = t(t(z)-zm)
      }

      Zc <- Re(apply(z, 2, FT))
      
      
      #Z <- Re(Zc)
      #Rc <- Re(apply(Zc, 2, function(col) wd_vec(col, p0, filter.number, family="DaubExPhase", bc="periodic")))
      Rc <- Re(apply(Zc, 2, wd_vec, p0, filter.number, family, bc))
      
      #R <- Re(Rc)
      m1=apply(abs(Zc),2,max)
      m2=apply(abs(Rc),2,max)
      temp=rbind(m1,m2)
      nr=apply(temp,2,max)
      ty <- rbind(Zc, Rc)
      nr = apply(abs(ty),2,max)
      quantiles_estimates = nr
    }
    
    if (type == "sqrt-LASSO") {
      quantiles_estimates <- numeric(m)
      
      z2 = matrix(rnorm(n*m,sd=sd), ncol=m)
      
      if (intercept){
        zm <- apply(z2,2,mean)
        z2  = t(t(z2)-zm)
      }

      
      Zc2 <- Re(apply(z2, 2, FT))
      #Z2 <- Re(Zc2)
      Rc2 <- Re(apply(Zc2, 2, function(col) wd_vec(col, p0, filter.number, family="DaubExPhase", bc="periodic")))
      #R2 <- Re(Rc2)
      ty2 <- rbind(Zc2, Rc2)
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
