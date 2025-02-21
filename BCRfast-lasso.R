# Définition de la fonction FISTA pour le Lasso
BCRfast_lasso = function( y, lambda, intercept=T, savecost=F, tol = 1e-20) {
  start <- Sys.time()
  p0=4
  filter.number=1
  J=log2(p0)
  
  p1 <- length(y)
  p <- 2*p1
  
  i0=as.numeric(intercept)
  beta <- rep(0,times=i0+p)
  
  continue=TRUE
  
  if(savecost) cost=lasso_cost(X, y, beta, lambda, type="LASSO", intercept=intercept)
  else cost=NA
  yt0 <- wd(rep(0,p1),filter.number=filter.number, family="DaubExPhase", bc="periodic")
  while(continue == TRUE){
    beta_old = beta
    
    u1 <- y - Re(FTinv(wr_vec(beta[-c((1:(i0+p1)))], p0=p0, yt0, filter.number=filter.number, family="DaubExPhase", bc="periodic")))
    if(intercept) u1=u1-beta[1]
    beta[i0+(1:p1)] <- soft_thresholding(Re(FT(u1)), lambda)

    u2 <- y - Re(FTinv(beta[i0+(1:p1)]))
    if(intercept) u2=u2-beta[1]
    beta[-c(1:(i0+p1))] = soft_thresholding(wd_vec(Re(FT(u2)), p0, filter.number, family="DaubExPhase", bc="periodic"), lambda)
    
    if(intercept){
      u3=y-Re(FTinv(wr_vec(beta[-c((1:(i0+p1)))], p0=p0, yt0, filter.number=filter.number, family="DaubExPhase", bc="periodic"))) - Re(FTinv(beta[i0+(1:p1)]))
      beta[1]=mean(u3)
    }
    
    # Vérification de la convergence)
    continue = ((sum((beta - beta_old)^2)/(1+sum(beta^2))) > tol) 
    # Mises à jour pour l'itération suivante
    if(savecost){
      newcost=lasso_cost(X, y, beta, lambda, type="LASSO", intercept=intercept)
      cost=c(cost,newcost)
    }
    #print(newcost)
     
  }
  end <- Sys.time()
  
  out=NULL
  out$beta=beta
  out$cost=cost
  out$time = end - start
  return(out)
}



