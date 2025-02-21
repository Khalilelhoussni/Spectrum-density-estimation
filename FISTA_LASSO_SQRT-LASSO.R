

fista_lasso_sqrt_lasso <- function(X, y, lambda, intercept=T,type ="LASSO" ,tol = 1e-20,rho=1.2,L0=0.05){
  # rho and L0 are linesearch parameters
  # m is monte-carlo parameter for lambdaqut
  startTime <- Sys.time()
  n <- nrow(X)
  p <- ncol(X)
  i0=as.numeric(intercept)
  beta <- rep(0,times=i0+p)
  
  u <- rep(0,times=p)
  t <- 1
  print("Calculation of L")
  L_lasso <- max(eigen(t(X)%*%X)$values)
  print("done")
  
  continue=T
  if (type =="LASSO"){

    cost=lasso_cost(X, y, beta, lambdaQUT3, type="LASSO", intercept=intercept)
    while(continue){
      
      beta_old = beta
      
      
      # Gradient step
      r=y-beta[1]*i0
      gradient <- - t(X) %*% (r - X %*% u)
      beta_temp <- u - (1 / L_lasso) * gradient
      
      # Proximal operator (soft-thresholding)
      beta[i0+(1:p)] <- soft_thresholding(beta_temp, lambdaQUT3/L_lasso)
      
      # Update t and z
      t_new <- (1 + sqrt(1 + 4 * t^2)) / 2
      u <- beta[i0+(1:p)] + ((t - 1) / t_new) * (beta[i0+(1:p)] - beta_old[i0+(1:p)])
      if (intercept){
        r <- y - X%*%beta[2:(p+1)]
        beta[1]<- mean(r)
      }
      
      # Vérification de la convergence)
      continue = ((sum((beta - beta_old)^2)/(1+sum(beta^2))) > tol) 
      # Mises à jour pour l'itération suivante
      cost=c(cost,lasso_cost(X, y, beta, lambdaQUT, type="LASSO", intercept=intercept))
      
      t <- t_new
    }
  }
  if (type =="sqrt-LASSO"){
    lambdaMC_sqrtLASSO=lambda_QUT(X,m,type="sqrt-LASSO",sd=sd,intercept = intercept)
    lambdaQUT3=lambdaMC_sqrtLASSO$lambdaQUT
    cost=lasso_cost(X, y, beta, lambdaQUT3, type="sqrt-LASSO", intercept=intercept)
    L <- L0
    while(continue){
      u1 <- u
      r=y-beta[1]*i0
      Xbeta=X %*% u
      gradient <- -t(X) %*% (r - Xbeta) / sqrt(sum((r - Xbeta)^2))
      beta_ls <- u - (1 / L) * gradient
      
      Pl <- rep(0,times=i0+p)
      Pl[i0+(1:p)] <- soft_thresholding(beta_ls, lambdaQUT3/L)
      if (intercept) {
        Pl[1] <- beta[1]
        u1 <- c(beta[1],u)
      }
      
      
      while (lasso_cost(X, y, Pl, lambdaQUT3, type="sqrt-LASSO", intercept=intercept) > Q_sqrtLASSO(X,y,Pl,u1,L,lambdaQUT3,intercept=intercept)){
        L <- rho*L
        Pl[i0+(1:p)] <- soft_thresholding(beta_ls, lambdaQUT3/L)
      }
      print(L)
      beta_old <- beta
      # Gradient step
      r=y-beta[1]*i0
      Xbeta=X %*% u
      gradient <- -t(X) %*% (r - Xbeta) / sqrt(sum((r - Xbeta)^2))
      beta_temp <- u - (1 / L) * gradient
      
      # Proximal operator (soft-thresholding)
      beta[i0+(1:p)] <- soft_thresholding(beta_temp, lambdaQUT3/L)
      
      # Update t and z
      t_new <- (1 + sqrt(1 + 4 * t^2)) / 2
      u <- beta[i0+(1:p)] + ((t - 1) / t_new) * (beta[i0+(1:p)] - beta_old[i0+(1:p)])
      if (intercept){
        r <- y - X%*%beta[2:(p+1)]
        beta[1]<- mean(r)
      }
      # Vérification de la convergence)
      continue = ((sum((beta - beta_old)^2)/(1+sum(beta^2))) > tol) 
      
      # Mises à jour pour l'itération suivante
      cost=c(cost,lasso_cost(X, y, beta, lambdaQUT3, type="sqrt-LASSO", intercept=intercept))
      print(lasso_cost(X, y, beta, lambdaQUT3, type="sqrt-LASSO", intercept=intercept))
      
      t <- t_new
    }
  }
  end <- Sys.time()
  out=NULL
  out$beta=beta
  out$cost = cost
  out$time = end - startTime
  return(out)
}  
  
