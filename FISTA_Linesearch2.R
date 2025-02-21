fista_lasso_sqrt_lasso_fast3 <- function(y,lambda, intercept=T, type ="LASSO", p0=4, filter.number=1, family="DaubExPhase", bc="periodic", tol = 1e-20,rho=1.1,L0=0.1){
  # rho and L0 are linesearch parameters
  # m is monte-carlo parameter for lambdaqut
  startTime <- Sys.time()
  J=log2(p0)
  n <- length(y)
  p <- 2*n
  i0=as.numeric(intercept)
  beta <- rep(0,times=i0+p)
  # Initialize matrices
  yt0 <- wd(rep(0,n), filter.number=filter.number, family=family, bc=bc)
  u <- rep(0,times=p)
  t <- 1
  continue=T
  A <- numeric(p)
  Pl <- rep(0,times=i0+p)
  if (type =="LASSO"){
    cost=lasso_cost_fast(y, beta, lambda, type="LASSO", intercept=intercept)
    L <- L0
    while(continue){
      u1 <- u
      r <- y-beta[1]*i0
      Z1 <- Re(FTinv(u[1:n]))
      R1 <- Re(FTinv(wr_vec(u[-c(1:n)], p0=p0, yt0, filter.number=filter.number, family=family, bc=bc)))
      Xu <- Z1 + R1
      P <- r - Xu
      Zc <- Re(FT(P))
      Rc <- Re(wd_vec(Zc, p0, filter.number, family=family, bc=bc))
      
      gradient <- - c(Zc,Rc)
      beta_ls <- u - (1 / L) * gradient
      Pl[i0+(1:p)] <- soft_thresholding(beta_ls, lambda/L)
      if (intercept) {
        Pl[1] <- beta[1]
        u1 <- c(beta[1],u)
      }
      
      while (lasso_cost_fast(y, Pl, lambda , type="LASSO", intercept=intercept) > Q_fista(y,Pl,beta,L,lambda ,intercept=intercept,type ="LASSO")){

        L <- rho*L
        Pl[i0+(1:p)] <- soft_thresholding(beta_ls, lambda /L)
      }
      
      beta_old = beta
      # Gradient step
      r=y-beta[1]*i0
      beta_temp <- u - (1 / L) * gradient
      
      # Proximal operator (soft-thresholding)
      beta[i0+(1:p)] <- soft_thresholding(beta_temp, lambda/L)
      
      # Update t and z
      t_new <- (1 + sqrt(1 + 4 * t^2)) / 2
      u <- beta[i0+(1:p)] + ((t - 1) / t_new) * (beta[i0+(1:p)] - beta_old[i0+(1:p)])
      
      if (intercept){
        A <- beta[2:(p+1)]
        Z1 <- Re(FTinv(A[1:n]))
        R1 <- Re(FTinv(wr_vec(A[-c(1:n)], p0=p0, yt0, filter.number=filter.number, family=family, bc=bc)))
        Xbeta2 <- Z1 + R1
        r <- y - Xbeta2
        beta[1] <- mean(r)
      }
      
      # Vérification de la convergence)
      continue = ((sum((beta - beta_old)^2)/(1+sum(beta^2))) > tol) 
      # Mises à jour pour l'itération suivante
      cost=c(cost, lasso_cost_fast(y, beta, lambda, type="LASSO", intercept=intercept))
      print(lasso_cost_fast(y, beta, lambda, type="LASSO", intercept=intercept))
      t <- t_new
    }
  }
  if (type =="sqrt-LASSO"){
    cost=lasso_cost_fast(y, beta, lambda , type="sqrt-LASSO", intercept=intercept)
    L <- L0
    L1 <- L
    while(continue){
      u1 <- u
      r=y-beta[1]*i0
      Zn <- Re(FTinv(u[1:n]))
      Rc <- Re(FTinv(wr_vec(u[-c(1:n)], p0=p0, yt0, filter.number=filter.number, family=family, bc=bc)))
      Xz <- Zn + Rc
      
      Z1 <- Re(FTinv(u[1:n]))
      R1 <- Re(FTinv(wr_vec(u[-c(1:n)], p0=p0, yt0, filter.number=filter.number, family=family, bc=bc)))
      Xu <- Z1 + R1
      P <- r - Xu
      Zi <- Re(FT(P))
      Ri <- Re(wd_vec(Zi, p0, filter.number, family=family, bc=bc))
      
      gradient <- - c(Zi,Ri) / sqrt(sum((r - Xz)^2))
      beta_ls <- u - (1 / L1) * gradient
      
      Pl[i0+(1:p)] <- soft_thresholding(beta_ls, lambda/L1)
      if (intercept) {
        Pl[1] <- beta[1]
        u1 <- c(beta[1],u)
      }
      
      while (lasso_cost_fast(y, Pl, lambda , type="sqrt-LASSO", intercept=intercept) < Q_fista(y,Pl,u1,L1,lambda,intercept=intercept, type ="sqrt-LASSO")){
        L1 <- L1/2
        Pl[i0+(1:p)] <- soft_thresholding(beta_ls, lambda/L1)
      }
      L <- min(2*L1,L0)
      L1 <- L
      beta_old <- beta
      
      beta_temp <- u - (1 / L) * gradient
      
      # Proximal operator (soft-thresholding)
      beta[i0+(1:p)] <- soft_thresholding(beta_temp, lambda /L)
      
      # Update t and z
      t_new <- (1 + sqrt(1 + 4 * t^2)) / 2
      u <- beta[i0+(1:p)] + ((t - 1) / t_new) * (beta[i0+(1:p)] - beta_old[i0+(1:p)])
      
      if (intercept){
        A <- beta[2:(p+1)]
        Z1 <- Re(FTinv(A[1:n]))
        R1 <- Re(FTinv(wr_vec(A[-c(1:n)], p0=p0, yt0, filter.number=filter.number, family=family, bc=bc)))
        Xbeta2 <- Z1 + R1
        r <- y - Xbeta2
        beta[1] <- mean(r)
      }
      # Vérification de la convergence)
      continue = ((sum((beta - beta_old)^2)/(1+sum(beta^2))) > tol) 
      
      # Mises à jour pour l'itération suivante
      cost=c(cost,lasso_cost_fast(y, beta, lambda , type="sqrt-LASSO", intercept=intercept))
      print(lasso_cost_fast(y, beta, lambda , type="sqrt-LASSO", intercept=intercept))
      
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
