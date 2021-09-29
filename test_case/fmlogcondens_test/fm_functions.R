Draw.boxplotsM <- function(Obj)
{
  p0hat <- Obj$p0hat.SP
  Sensitivity <- Obj$Sensitivity
  FPR <- 1 - Obj$Specificity
  runtime <- Obj$t
  
  par(mfrow = c(1, 4))
  boxplot(p0hat, main = "Estimates of p0")
  boxplot(FPR, main = "False Positive Rate = (1 - Specificity)", ylim = c(0, 0.05))
  boxplot(Sensitivity, main = "Sensitivity", ylim = c(0.8, 1))
  boxplot(runtimetime, main="time")
}


SimMultNormal <- function(M, n, p0)
  # [Multivariate] Null: Normal, Nonnull: Normal
{
  library(copula)
  param <- 0.5
  
  library(mvtnorm)
  rho12 <- 0.25
  Sigma0 <- matrix(c(1, rho12, 
                     rho12, 1), 
                   ncol = 2, byrow = T)
  
  result <- data.frame(p0hat = rep(NA, M),
                       Sensitivity = rep(NA, M),
                       Specificity = rep(NA, M),
                       t=rep(NA,M))
  
  for ( r in 1:M ) {
    cat(r, "/", M, "\n")
    n0 <- rbinom(1, n, p0)
    n1 <- n - n0
    V <- rCopula(n1, ellipCopula(family = "normal", param = param, dim = 2))
    z0 <- rmvnorm(n0, sigma = Sigma0)
    z1 <- qnorm(V, mean = 3.5, sd = .5)
    z <- rbind(z0, z1)
    plot(z, col = c(rep(1, n0), rep(2, n1)), pch = 20)
    
    runtime<-system.time(res <- sp.mix.multi(z))[3]
    
    p0hat <- res$p.0
    Nhat <- as.integer(res$localfdr <= 0.2)
    TP <- sum(Nhat[-(1:n0)] == 1)
    TN <- sum(Nhat[1:n0] == 0)
    FP <- n0 - TN
    FN <- n1 - TP
    result$p0hat.SP[r] <- p0hat
    result$Sensitivity[r] <- TP/(TP + FN)
    result$Specificity[r] <- TN/(TN + FP)
    result$t[r] <- runtime
  }
  
  return(result)
}

SimMultNormal3d <- function(M, n, p0)
  # [Multivariate] Null: Normal, Nonnull: Normal
{
  library(copula)
  param <- 0.5
  
  library(mvtnorm)
  rho12 <- 0.25
  rho13 <- 0.25
  rho23 <- 0.25
  Sigma0 <- matrix(c(1, rho12, rho13,
                     rho12, 1, rho23,
                     rho13,rho23,1), 
                   ncol = 3, byrow = T)
  
  result <- data.frame(p0hat = rep(NA, M),
                       Sensitivity = rep(NA, M),
                       Specificity = rep(NA, M),
                       t=rep(NA,M))
  
  for ( r in 1:M ) {
    cat(r, "/", M, "\n")
    n0 <- rbinom(1, n, p0)
    n1 <- n - n0
    V <- rCopula(n1, ellipCopula(family = "normal", param = param, dim = 3))
    z0 <- rmvnorm(n0, sigma = Sigma0)
    z1 <- qnorm(V, mean = 3.5, sd = .5)
    z <- rbind(z0, z1)
    
    runtime<-system.time(res <- sp.mix.multi(z))[3]
    
    p0hat <- res$p.0
    Nhat <- as.integer(res$localfdr <= 0.2)
    TP <- sum(Nhat[-(1:n0)] == 1)
    TN <- sum(Nhat[1:n0] == 0)
    FP <- n0 - TN
    FN <- n1 - TP
    result$p0hat.SP[r] <- p0hat
    result$Sensitivity[r] <- TP/(TP + FN)
    result$Specificity[r] <- TN/(TN + FP)
    result$t[r] <- runtime
  }
  
  return(result)
}


SimMultGamma <- function(M, n, p0)
  # [Multivariate] Null: Normal, Nonnull: Normal
{
  library(copula)
  param <- 2
  
  rho12 <- 0.5
  Sigma0 <- matrix(c(1, rho12, 
                     rho12, 1), 
                   ncol = 2, byrow = T)
  
  library(mvtnorm)
  
  result <- data.frame(p0hat.SP = rep(NA, M),
                       Sensitivity = rep(NA, M),
                       Specificity = rep(NA, M),
                       t=rep(NA,M))
  
  for ( r in 1:M ) {
    cat(r, "/", M, "\n")
    n0 <- rbinom(1, n, p0)
    n1 <- n - n0
    V <- rCopula(n1, frankCopula(param = param, dim = 2))
    z0 <- rmvnorm(n0, sigma = Sigma0)
    z1 <- qgamma(V, shape = 12, rate = 4)
    z <- rbind(z0, z1)
    plot(z, col = c(rep(1, n0), rep(2, n1)), pch = 20)
    
    runtime<-system.time(res <- sp.mix.multi(z))[3]
    p0hat <- res$p.0
    Nhat <- as.integer(res$localfdr <= 0.2)
    TP <- sum(Nhat[-(1:n0)] == 1)
    TN <- sum(Nhat[1:n0] == 0)
    FP <- n0 - TN
    FN <- n1 - TP
    result$p0hat.SP[r] <- p0hat
    result$Sensitivity[r] <- TP/(TP + FN)
    result$Specificity[r] <- TN/(TN + FP)
    result$t[r] <- runtime
    
  }
  
  return(result)
}

SimMultGamma3d <- function(M, n, p0)
  # [Multivariate] Null: Normal, Nonnull: Normal
{
  library(copula)
  param <- 2
  
  rho12 <- 0.25
  rho13 <- 0.25
  rho23 <- 0.25
  Sigma0 <- matrix(c(1, rho12, rho13,
                     rho12, 1, rho23,
                     rho13,rho23,1), 
                   ncol = 3, byrow = T)
  
  library(mvtnorm)
  
  result <- data.frame(p0hat.SP = rep(NA, M),
                       Sensitivity = rep(NA, M),
                       Specificity = rep(NA, M),
                       t=rep(NA,M))
  
  
  for ( r in 1:M ) {
    cat(r, "/", M, "\n")
    n0 <- rbinom(1, n, p0)
    n1 <- n - n0
    V <- rCopula(n1, frankCopula(param = param, dim = 3))
    z0 <- rmvnorm(n0, sigma = Sigma0)
    z1 <- qgamma(V, shape = 12, rate = 4)
    z <- rbind(z0, z1)
    runtime<-system.time(res <- sp.mix.multi(z))[3]
    p0hat <- res$p.0
    Nhat <- as.integer(res$localfdr <= 0.2)
    TP <- sum(Nhat[-(1:n0)] == 1)
    TN <- sum(Nhat[1:n0] == 0)
    FP <- n0 - TN
    FN <- n1 - TP
    result$p0hat.SP[r] <- p0hat
    result$Sensitivity[r] <- TP/(TP + FN)
    result$Specificity[r] <- TN/(TN + FP)
    result$t[r] <- runtime
    
  }
  
  return(result)
}

