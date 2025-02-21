# BCR to solve LASSO with or without an intercept
BCR_lasso = function(X1, X2, y, lambda, intercept=T, tol = 1e-20) {
  start <- Sys.time()
  
  X = cbind(X1,X2)
  n <- nrow(X)
  p <- ncol(X)
  p1 <- ncol(X1)
  i0=as.numeric(intercept)
  beta <- rep(0,times=i0+p)

  continue=TRUE
  cost=lasso_cost(X, y, beta, lambda, type="LASSO", intercept=intercept)

  while(continue == TRUE){
    beta_old = beta

    u1 <- y - X2%*%beta[-c((1:(i0+p1)))]
    if(intercept) u1=u1-beta[1]
    beta[i0+(1:p1)] <- soft_thresholding(t(X1)%*% u1, lambda)
    
    u2 <- y - X1%*%beta[i0+(1:p1)]
    if(intercept) u2=u2-beta[1]
    beta[-c(1:(i0+p1))] = soft_thresholding(t(X2) %*% u2, lambda)
    if(intercept){
      u3=y-X%*%beta[-1]
      beta[1]=mean(u3)
    }
    
    continue = ((sum((beta - beta_old)^2)/(1+sum(beta^2))) > tol) 
    # Vérification de la convergence)
    # Mises à jour pour l'itération suivante
    newcost=lasso_cost(X, y, beta, lambda, type="LASSO", intercept=intercept)
    cost=c(cost,newcost)

  }
  end <- Sys.time()
  
  out=NULL
  out$beta=beta
  out$cost=cost
  out$time = end - start
  return(out)
}





