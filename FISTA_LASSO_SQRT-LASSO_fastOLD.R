fista_lasso_sqrt_lasso_fast <- function( y, lambda, intercept=T, type ="LASSO", p0=4, filter.number=1, family="DaubExPhase", bc="periodic", tol = 1e-20,rho=1.2,L0=0.05){
  # rho and L0 are linesearch parameters
  # m is monte-carlo parameter for lambdaqut
  startTime <- Sys.time()
  J=log2(p0)
  n <- length(y)
  p <- 2*n
  i0=as.numeric(intercept)
  beta <- rep(0,times=i0+p)
  
  Zc <- FT(y)
  Z  <- Re(Zc)
  Rc <- wd_vec(Z, p0, filter.number, family=family, bc=bc)
  R  <- Re(Rc)
  Xy <- c(Z,R)
  
  
  
  yt0=wd(rep(0,n), filter.number=filter.number, family=family, bc=bc)
  T2 = matrix(NA, n, n)
  for(i in 1:n){
    yy=wr_vec(diag(n)[,i], p0=p0, yt0, filter.number=filter.number, family="DaubExPhase", bc="periodic")
    T2[,i]=yy
  }
  
  T1 =matrix(NA, n, n)
  for(i in 1:n){
    h<- diag(n)[,i]
    yX1=wd_vec(h, p0, filter.number, family="DaubExPhase", bc="periodic")
    T1[,i]=yX1
  }
  
  
  
  a1 <- cbind(diag(n),T2)
  a2 <- cbind(T1,diag(n))
  Xx <- rbind(a1,a2)
  
  u <- rep(0,times=p)
  t <- 1
  print("Calculation of L")
  L_lasso <- max(eigen(Xx)$values)
  print("done")
  
  continue=T
  if (type =="LASSO"){
    cost=lasso_cost_fast(y, beta, lambda, type="LASSO", intercept=intercept)
    k<- numeric(p)
    while(continue){
      
      beta_old = beta
      # Gradient step
      r=y-beta[1]*i0
      Zc <- FT(r)
      Z  <- Re(Zc)
      Rc <- wd_vec(Z, p0, filter.number, family=family, bc=bc)
      R  <- Re(Rc)
      Xr <- c(Z,R)
      gradient <- (-Xr + Xx %*% u)
      beta_temp <- u - (1 / L_lasso) * gradient
      
      # Proximal operator (soft-thresholding)
      beta[i0+(1:p)] <- soft_thresholding(beta_temp, lambda/L_lasso)
      
      # Update t and z
      t_new <- (1 + sqrt(1 + 4 * t^2)) / 2
      u <- beta[i0+(1:p)] + ((t - 1) / t_new) * (beta[i0+(1:p)] - beta_old[i0+(1:p)])
      
      if (intercept){
        k <- beta[2:(p+1)]
        Z1 <- FTinv(k[1:n])
        Z2 <- Re(Z1)
        R1 <- FTinv(wr_vec(k[-c(1:n)], p0=p0, yt0, filter.number=filter.number, family=family, bc=bc))
        R2 <- Re(R1)
        Xbeta2 <- Z2 + R2
        r <- y - Xbeta2
        beta[1] <- mean(r)
      }
      
      # Vérification de la convergence)
      continue = ((sum((beta - beta_old)^2)/(1+sum(beta^2))) > tol) 
      # Mises à jour pour l'itération suivante
      cost=c(cost, lasso_cost_fast( y, beta, lambda, type="LASSO", intercept=intercept))
      
      t <- t_new
    }
  }
  if (type =="sqrt-LASSO"){
    cost=lasso_cost_fast(y, beta, lambda , type="sqrt-LASSO", intercept=intercept)
    L <- L0
    while(continue){
      u1 <- u
      r=y-beta[1]*i0
      Zn <- FTinv(u[1:n])
      Zm <- Re(Zn)
      Rc <- FTinv(wr_vec(u[-c(1:n)], p0=p0, yt0, filter.number=filter.number, family=family, bc=bc))
      R <- Re(Rc)
      Xz <- Zm + R
      
      Zc <- Re(FT(r))
      #Z  <- Re(Zc)
      Rc2 <- Re(wd_vec(Z, p0, filter.number, family=family, bc=bc))
      #R2  <- Re(Rc2)
      Xr <- c(Zc,Rc2)
      
      gradient <- (-Xr + Xx %*% u) / sqrt(sum((r - Xz)^2))
      beta_ls <- u - (1 / L) * gradient
      
      Pl <- rep(0,times=i0+p)
      Pl[i0+(1:p)] <- soft_thresholding(beta_ls, lambda/L)
      if (intercept) {
        Pl[1] <- beta[1]
        u1 <- c(beta[1],u)
      }
      
      while (lasso_cost_fast(y, Pl, lambda , type="sqrt-LASSO", intercept=intercept) > Q_fista(y,Pl,u1,L,lambda ,intercept=intercept, type = "sqrt-LASSO")){
        print(lasso_cost_fast(y, Pl, lambda , type="sqrt-LASSO", intercept=intercept))
        print(Q_fista(y,Pl,u1,L,lambda ,intercept=intercept))
        L <- rho*L
        Pl[i0+(1:p)] <- soft_thresholding(beta_ls, lambda /L)
        
      }
      
      beta_old <- beta
      # Gradient step
      r=y-beta[1]*i0
      gradient <- (-Xr + Xx %*% u) / sqrt(sum((r - Xz)^2))
      beta_temp <- u - (1 / L) * gradient
      
      # Proximal operator (soft-thresholding)
      beta[i0+(1:p)] <- soft_thresholding(beta_temp, lambda /L)
      
      # Update t and z
      t_new <- (1 + sqrt(1 + 4 * t^2)) / 2
      u <- beta[i0+(1:p)] + ((t - 1) / t_new) * (beta[i0+(1:p)] - beta_old[i0+(1:p)])
      if (intercept){
        k <- beta[2:(p+1)]
        Z1 <- FTinv(k[1:n])
        Z2 <- Re(Z1)
        R1 <- FTinv(wr_vec(k[-c(1:n)], p0=p0, yt0, filter.number=filter.number, family=family, bc=bc))
        R2 <- Re(Rc)
        Xbeta2 <- Z2 + R2
        r <- y - Xbeta2
        beta[1] <- mean(r)
      }
      # Vérification de la convergence)
      continue = ((sum((beta - beta_old)^2)/(1+sum(beta^2))) > tol) 
      
      # Mises à jour pour l'itération suivante
      cost=c(cost,lasso_cost_fast( y, beta, lambda , type="sqrt-LASSO", intercept=intercept))
      print(lasso_cost_fast( y, beta, lambda , type="sqrt-LASSO", intercept=intercept))
      
      t <- t_new
    }
  }
  end <- Sys.time()
  out=NULL
  out$beta=beta
  out$cost = cost
  out$time = end - startTime
  out$lambda = lambda
  return(out)
}
