fista_nonconvex <- function(X, y, lambda, intercept=T, tol = 1e-8,a=1) {
  startTime <- Sys.time()
  n <- nrow(X)
  p <- ncol(X)
  i0=as.numeric(intercept)
  beta <- rep(0,times=i0+p)
  
  z <- rep(0,times=p)
  t <- 1
  print("Calculation of L")
  L <- max(eigen(t(X)%*%X)$values)
  print("done")
  cost=lasso_cost_NC2(X, y, beta, lambda, intercept=intercept)
  
  continue=T
  
  
  while(continue){
    
    beta_old = beta
    
    
    # Gradient step
    r=y-beta[1]*i0
    gradient <- - t(X) %*% (r - X %*% z)
    beta_temp <- z - (1 / L) * gradient
    
    # Proximal operator (soft-thresholding)
    beta[i0+(1:p)] <- NC_Thresholding3(beta_temp, lambda/L,a)
    
    # Update t and z
    t_new <- (1 + sqrt(1 + 4 * t^2)) / 2
    z <- beta[i0+(1:p)] + ((t - 1) / t_new) * (beta[i0+(1:p)] - beta_old[i0+(1:p)])
    if (intercept){
      r <- y - X%*%beta[2:(p+1)]
      beta[1]<- mean(r)
    }
    
    # Vérification de la convergence)
    continue = ((sum((beta - beta_old)^2)/(1+sum(beta^2))) > tol) 
    # Mises à jour pour l'itération suivante
    cost=c(cost,lasso_cost_NC2(X, y, beta, lambda, intercept=intercept))
    
    t <- t_new
  }
  end <- Sys.time()
  out=NULL
  out$beta=beta
  out$cost = cost
  out$time = end - startTime
  return(out)
}

