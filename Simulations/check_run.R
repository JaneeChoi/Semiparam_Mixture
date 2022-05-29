library(SpMix)
set.seed(990605)

M=1
n=1000
p0=0.8

SimMultNormal <- function(M, n, p0)
  # [Multivariate] Null: Normal, Nonnull: Normal
{
  library(copula)
  param <- 0.5
  
  library(mclust)
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
    result$p0hat[r] <- p0hat
    result$Sensitivity[r] <- TP/(TP + FN)
    result$Specificity[r] <- TN/(TN + FP)
    result$t[r] <- runtime
  }
  return(result)
}


res<-SimMultNormal(1,1000,0.8)
print(res)
