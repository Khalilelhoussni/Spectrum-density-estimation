


r <- function(x, lambda, a) { 
  r2_val <- 0.5 * (x - a) + sqrt(abs(0.25 * (x + a)^2 - lambda))
  print(log(1 + abs(r2_val)/a))
  r3 <- 0.5 * (x - r2_val)^2 + lambda * log(1 + abs(r2_val)/a) - x^2/2
  return(r3)
}

NC_Thresholding3 <- function(x, lambda, a) {
  if (sqrt(lambda) <= a) {
    prox <- numeric(length(x))
    prox[abs(x) <= lambda/a] <- 0
    idx <- abs(x) > lambda/a
    r2 <- 0.5 * (abs(x[idx]) - a) + sqrt(0.25 * (abs(x[idx]) + a)^2 - lambda)
    prox[idx] <- sign(x[idx]) * r2
  } else {
    prox <- numeric(length(x))
    prox[abs(x) < 2 * sqrt(lambda) - a] <- 0
    idx <- abs(x) >= 2 * sqrt(lambda) - a

    
    Phi0 <- 2*sqrt(lambda) - a
    Phi1 <-  lambda/a
    x1 <- uniroot(r, c(Phi0, Phi1), lambda = lambda, a = a)$root
    prox[idx] <- ifelse(abs(x[idx]) <= x1, 0, sign(x[idx]) * (0.5 * (abs(x[idx]) - a) + sqrt(0.25 * (abs(x[idx]) + a)^2 - lambda)))
  }
  return(prox)
}

