#source(file="SPMix-LocalFDR/SpMix.R")
#source(file = 'test_case/fmlogcondens_test/sp.mix.multi.fmlogcondens.R')

devtools::install_github("JaneeChoi/SpMix")
library(SpMix)
set.seed(210828)

Draw.boxplotsM <- function(Obj)
{
  p0hat <- Obj$p0hat.SP
  Sensitivity <- Obj$Sensitivity
  FPR <- 1 - Obj$Specificity
  
  par(mfrow = c(1, 3))
  boxplot(p0hat, main = "Estimates of p0")
  boxplot(FPR, main = "False Positive Rate = (1 - Specificity)", ylim = c(0, 0.05))
  boxplot(Sensitivity, main = "Sensitivity", ylim = c(0.8, 1))
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
                       Specificity = rep(NA, M))
  
  for ( r in 1:M ) {
    cat(r, "/", M, "\n")
    n0 <- rbinom(1, n, p0)
    n1 <- n - n0
    V <- rCopula(n1, ellipCopula(family = "normal", 
                                 param = param, dim = 2))
    z0 <- rmvnorm(n0, sigma = Sigma0)
    z1 <- qnorm(V, mean = 3.5, sd = .5)
    z <- rbind(z0, z1)
    plot(z, col = c(rep(1, n0), rep(2, n1)), pch = 20)

    res <- sp.mix.multi(z)
    
    p0hat <- res$p.0
    Nhat <- as.integer(res$localfdr <= 0.2)
    TP <- sum(Nhat[-(1:n0)] == 1)
    TN <- sum(Nhat[1:n0] == 0)
    FP <- n0 - TN
    FN <- n1 - TP
    result$p0hat.SP[r] <- p0hat
    result$Sensitivity[r] <- TP/(TP + FN)
    result$Specificity[r] <- TN/(TN + FP)
  }
  
  return(result)
}


M <- 1


for (N in c(100,500,1000,5000,10000,15000,20000,25000)){
  time1<-system.time(Res.1<-SimMultNormal(M = M, n = N, p0 = 0.95))[3]
  print("finish 1")
  time2<-system.time(Res.2<-SimMultNormal(M = M, n = N, p0 = 0.90))[3]
  print("finish 2")
  time3<-system.time(Res.3<-SimMultNormal(M = M, n = N, p0 = 0.80))[3]
  print("finish 3")
  result_df<-data.frame(res1<-c(time1,Res.1$p0hat.SP[1],Res.1$Sensitivity[1],Res.1$Specificity[1]),res2<-c(time2,Res.2$p0hat.SP[1],Res.2$Sensitivity[1],Res.2$Specificity[1]),res3<-c(time3,Res.3$p0hat.SP[1],Res.3$Sensitivity[1],Res.3$Specificity[1]))
  write.csv(result_df,paste0("result_",N,".csv"))
  print("finish")
  print(N)
}

# sp.mix.multi using fmlogcondens

devtools::install_github("JaneeChoi/SpMix",ref = "fmlogcondens")
library(SpMix)
set.seed(210828)

M <- 1
setwd("/Users/choiiiiii/Documents/GitHub/Semiparam_Mixture/test_case/fmlogcondens_test/Multi_Normal_fmlogcondens")
for (N in c(100,500,1000,5000,10000,15000,20000,25000)){
  time1<-system.time(Res.1<-SimMultNormal(M = M, n = N, p0 = 0.95))[3]
  print("finish 1")
  time2<-system.time(Res.2<-SimMultNormal(M = M, n = N, p0 = 0.90))[3]
  print("finish 2")
  time3<-system.time(Res.3<-SimMultNormal(M = M, n = N, p0 = 0.80))[3]
  print("finish 3")
  result_df<-data.frame(res1<-c(time1,Res.1$p0hat.SP[1],Res.1$Sensitivity[1],Res.1$Specificity[1]),res2<-c(time2,Res.2$p0hat.SP[1],Res.2$Sensitivity[1],Res.2$Specificity[1]),res3<-c(time3,Res.3$p0hat.SP[1],Res.3$Sensitivity[1],Res.3$Specificity[1]))
  write.csv(result_df,paste0("result_",N,".csv"))
  print("finish")
  print(N)
}

# sp.mix.multi using fmlogcondens

devtools::install_github("JaneeChoi/SpMix",ref = "fmlogcondens")
library(SpMix)
set.seed(210828)

M <- 1
setwd("/Users/choiiiiii/Documents/GitHub/Semiparam_Mixture/test_case/fmlogcondens_test/Multi_Gamma_fmlogcondens")
for (N in c(100,500,1000,5000,10000,15000,20000,25000)){
  time1<-system.time(Res.1<-SimMultNormal(M = M, n = N, p0 = 0.95))[3]
  print("finish 1")
  time2<-system.time(Res.2<-SimMultNormal(M = M, n = N, p0 = 0.90))[3]
  print("finish 2")
  time3<-system.time(Res.3<-SimMultNormal(M = M, n = N, p0 = 0.80))[3]
  print("finish 3")
  result_df<-data.frame(res1<-c(time1,Res.1$p0hat.SP[1],Res.1$Sensitivity[1],Res.1$Specificity[1]),res2<-c(time2,Res.2$p0hat.SP[1],Res.2$Sensitivity[1],Res.2$Specificity[1]),res3<-c(time3,Res.3$p0hat.SP[1],Res.3$Sensitivity[1],Res.3$Specificity[1]))
  write.csv(result_df,paste0("result_gamma_",N,".csv"))
  print("finish")
  print(N)
}



# Draw boxplot

png(file = "MultiNormalp95.png", height = 600, width = 900)
Draw.boxplotsM(Res.7)
dev.off()

#png(file = "MultiNormalp90.png", height = 600, width = 900)
Draw.boxplotsM(Res.8)
#dev.off()

#png(file = "MultiNormalp80.png", height = 600, width = 900)
Draw.boxplotsM(Res.9)
#dev.off()

