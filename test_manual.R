## First, we need package 'devtools', 'LogConcDEAD'

install.packages(c("devtools","LogConcDEAD"))

## Install JaneeChoi/fmlogcondens

devtools::install_github("JaneeChoi/fmlogcondens")

# simple check for JaneeChoi/fmlogcondens

library(fmlogcondens)

# draw samples from normal distribution
X <- matrix(rnorm(500),250,2)
# estimate log-concave density
system.time(r <- fmlcd(X))

# load `LogConcDEAD` library for plotting capatibilities
library(LogConcDEAD)
r <- LogConcDEAD::getinfolcd(X, r$logMLE) # create a `LogConcDEAD` object
# plot estimated density
par(mfrow = c(1, 2)) #square plots
plot(r, addp = FALSE, asp = 1, main="density")
plot(r, uselog = TRUE, addp = FALSE, asp = 1, main="log density")

# draw samples from normal distribution
X <- matrix(rnorm(200),100,2) 
# calculate parameters of convex hull of X
r <- calcCvxHullFaces(X)
# draw random parameters of 10 hyperplanes
a <- matrix(runif(10*2),10,2)
b <- runif(10)

# calculate parameters of convex hull of X
correctIntegral(X,rep(0,2),a,b,r$cvh)

# find initial hyperplane parameters based on a smooth log-concave density
paramFitGammaOne(X, rep(1 / nrow(X), nrow(X)), r$ACVH, r$bCVH, r$cvh)

# find initial hyperplane parameters based on a kernel density estimator with Gaussian kernel
paramFitKernelDensity(X, rep(1 / nrow(X), nrow(X)), r$cvh)

## Install JaneeChoi/SpMix which imports JaneeChoi/fmlogcondens

devtools::install_github("JaneeChoi/SpMix",force=TRUE)

# simple runtime check for JaneeChoi/SpMix

library(SpMix)

set.seed(222)
X = matrix(rnorm(2000), 1000, 2)

system.time(r <- fmlcd(X))
str(r) 

# runtime check for 2d normal mixture (M=1)

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
    result$p0hat[r] <- p0hat
    result$Sensitivity[r] <- TP/(TP + FN)
    result$Specificity[r] <- TN/(TN + FP)
    result$t[r] <- runtime
  }
  
  return(result)
}

result_df=data.frame(N=c(),p0=c(),time=c())

for (N in c(100,500,1000,5000,10000,15000,20000,25000)){
  for (p in c(0.95,0.90,0.80)){
    Res<-SimMultNormal(M = 1, n = N, p0 = p)
    result_df=rbind(result_df,data.frame(N=N,p0=p,time=Res$t))
  }
}

# plot and save the results

library(ggplot2)
library(dplyr)

p1<-ggplot(data=result_df,aes(x=N,y=time,group=p0,colour=p0))+geom_line()
ggsave("time.png", plot=p1, height=6, width=8, dpi=600)


