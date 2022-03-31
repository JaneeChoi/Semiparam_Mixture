rv<-read.table("/Users/user/Documents/GitHub/Semiparam_Mixture/Galaxy/rvdata_for_stat.txt", skip = 1)
z_vel = scale(rv$V2)[,1]
length(z_vel)
hist(z_vel,breaks=24,main="z-value for velocity",xlab="z-value")

dim(z)[2]
library(multiLocalFDR)
library(LogConcDEAD)


library(mclust)
library(logcondens)

SpMix <- function(z, tol = 5e-6, alternative = "greater", max_iter = 30, mono = TRUE, thre_z = 0.9,
                  Uthre_gam = 0.9, Lthre_gam = 0.01 )
{
  # *****************DEFINITION OF INTERNAL FUNCTIONS ******************

  NormMix <- function(z, tol = 5e-3, max_iter = 10)
  {

    k <- 0; converged <- 0
    z <- as.matrix(z)
    m_dist <- mahalanobis(z, 0, cov(z))

    p0 <- mean(m_dist <= 1.65)
    mu0 <- rep(0, dim(z)[2])
    sig0 <- diag(1, dim(z)[2])
    f0 <- dmvnorm(z, mu0, sig0)
    mu1 <- apply(as.matrix(z[m_dist > 1.65,]), 2, mean)
    sig1 <- cov(as.matrix(z[m_dist > 1.65,]))
    f1 <- dmvnorm(z, mu1, sig1)

    while ((k < 3) | ((k < max_iter) & (!converged))) {
      k <- k + 1

      ## E-step
      gam <- p0 * f0 / (p0 * f0 + (1-p0) * f1)

      ## M-step
      new_p0 <- mean(gam)
      new_mu0 <- as.vector(t(z) %*% gam) / sum(gam)
      dev0 <- (z - new_mu0) * sqrt(gam)
      new_sig0 <- t(dev0) %*% dev0 / sum(gam)
      f0 <- dmvnorm(z, new_mu0, new_sig0)
      new_mu1 <- as.vector(t(z) %*% (1 - gam)) / sum(1 - gam)
      dev1 <- (z - new_mu1) * sqrt(1 - gam)
      new_sig1 <- t(dev1) %*% dev1 / sum(1 - gam)

      ## Update
      diff <- max(abs((new_mu0 - mu0)),
                  abs((new_sig0 - sig0)),
                  abs((new_mu1 - mu1)),
                  abs((new_sig1 - sig1)),
                  abs((new_p0 - p0)))
      converged <- (diff <= tol)
      p0 <- new_p0
      mu0 <- new_mu0
      sig0 <- new_sig0
      mu1 <- new_mu1
      sig1 <- new_sig1
    }

    return(list(p0 = p0, mu0 = mu0, sig0 = sig0, mu1 = mu1, sig1 = sig1))
  }



  NE <- function(x, X)
  {
    n <- nrow(X)
    p <- ncol(X)
    xx <- matrix(x, nrow = n, ncol = p, byrow = TRUE)
    ne.ind <- apply(1*(X >= xx), 1, prod)

    return((1:n)[ne.ind == 1])
  }


  MonotoneFDR <- function(z, fdr)
  {
    n <- nrow(z)
    MFDR <- numeric(n)
    for (i in 1:n) {
      MFDR[i] <- max(fdr[NE(z[i,], z)])
    }

    return(MFDR)
  }

  # ******************* MAIN FUNCTION *******************************

  z <- as.matrix(z_vel)
  n <- dim(z)[1]
  print(n)

  ## Initial step: to fit normal mixture
  if (dim(z)[2] == 1) {
    if (alternative == "greater" | alternative == "g") {
      q0 <- quantile(z, probs = .9)
      p0 <- mean(z <= q0)
      mu0 <- mean(z[z <= q0])
      sig0 <- sd(z[z <= q0])
      f0 <- dmvnorm(z, mu0, sig0)
      mu1 <- mean(z[z > q0])
      sig1 <- sd(z[z > q0])
      f1 <- dnorm(z, mu1, sig1)
    }
    else {
      z_vec = z[,1]
      q0 <- quantile(z_vec, probs = .7)
      p0 <- mean(z_vec >= q0)
      mu0 <- mean(z_vec[z_vec >= q0])
      sig0 <- sd(z_vec[z_vec >= q0])
      f0 <- dmvnorm(z_vec, mu0, sig0)
      mu1 <- mean(z_vec[z_vec < q0])
      sig1 <- sd(z_vec[z_vec < q0])
      f1 <- dnorm(z_vec, mu1, sig1)
    }
  }
  else {
    Params <- NormMix(z)
    p0 <- Params$p0
    mu0 <- Params$mu0
    sig0 <- Params$sig0
    f1 <- dmvnorm(z, Params$mu1, Params$sig1)
  }
  gam <- f <- rep(0, n)

  ## EM-step
  k <- 0; converged <- 0
  while ( (k < 3)|((k < max_iter) & (!converged)) ) {
    k <- k + 1

    ## E-step
    new_f <- (p0 * f0 + (1 - p0) * f1)
    new_gam <- p0 * f0 / new_f

    ## M-step
    sum_gam <- sum(new_gam)
    new_mu0 <- as.vector(t(z) %*% new_gam) / sum_gam
    dev <- t(t(z)-new_mu0) * sqrt(new_gam)
    new_sig0 <- t(dev) %*% dev / sum_gam
    new_p0 <- mean(new_gam)
    new_f0 <- dmvnorm(z, new_mu0, new_sig0)
    weight <- 1 - new_gam
    new_f1 <- rep(0, n)
    which_z <- (new_gam <= thre_z)
    lcd <- mlelcd(x = z[,1][which_z], w = weight[which_z])
    new_f1[which_z]  <- exp(lcd$logMLE)[rank(z[which_z,])]

    ## Update
    which_gam <- (new_gam <= Uthre_gam) * (new_gam >= Lthre_gam)
    diff <- max(abs(gam - new_gam)[which_gam])
    converged <- (diff <= tol)
    cat("   EM iteration:", k, ", Change in mdfdr fit = ", round(diff, 5), "\n")
    p0 <- new_p0; mu_0 <- new_mu0; sig0 <- new_sig0
    f1 <- new_f1
    f0 <- new_f0
    f <- new_f
    gam <- new_gam
  }

  res <- list(p0 = p0, mu0 = mu0, sig0 = sig0,
              f = f, f1 = f1, localFDR = gam, iter = k)

  return(res)
}

SpMix_vel <- SpMix(z_vel, tol = 5e-6, alternative = "less",thre_z = 0.9,
                   Uthre_gam = 0.5, Lthre_gam = 0.01)


plotFDR(z_vel, SpMix_vel$p0, SpMix_vel$mu0, SpMix_vel$sig0, SpMix_vel$f1, SpMix_vel$localFDR, alternative = "less")


sp.mix.1D <- function(z, tol = 5.0e-6, max.iter = 30, doplot = TRUE, thre.localFDR = 0.2, thre.z = 0.95, Uthre.gam = 0.9, Lthre.gam = 0.01)
{
  #library(LogConcDEAD)
  library(logcondens)

  z <- as.numeric(z)
  n <- length(z)

  ## Initial step
  q0 <- quantile(z, probs = .9)
  p.0 <- mean(z >= q0)
  mu.0 <- mean(z[z >= q0])
  sig.0 <- sd(z[z >= q0])

  mu.1 <- mean(z[z < q0])
  sig.1 <- sd(z[z < q0])
  f1.tilde <- dnorm(z, mean = mu.1, sd = sig.1)

  f <- gam <- rep(0, n)

  ## EM-step
  k <- 0; converged <- 0
  while ( (k < 3) | ((k < max.iter) & (!converged)) ) {
    k <- k + 1
    ## E-step
    tmp <- p.0*dnorm(z, mu.0, sig.0)
    new.f <- tmp + (1 - p.0)*f1.tilde
    new.gam <- tmp/new.f

    ## M-step
    w.gam <- new.gam/sum(new.gam, na.rm = TRUE)
    new.mu.0 <- sum(w.gam*z, na.rm = TRUE)
    new.sig.0 <- sqrt(sum(w.gam*(z-new.mu.0)^2, na.rm = TRUE))
    new.p.0 <- mean(new.gam, na.rm = TRUE)

    new.f1.tilde <- rep(0, n)
    which.z <- new.gam <= thre.z
    weight <- 1 - new.gam[which.z]
    weight <- weight/sum(weight)
    lcd <- mlelcd(x = z[which.z], w = weight)
    new.f1.tilde[which.z] <- exp(lcd$logMLE)[rank(z[which.z])]
    which.gam <- (new.gam <= Uthre.gam)*(new.gam >= Lthre.gam)
    diff <- max(abs(gam - new.gam)[which.gam])
    converged <- (diff <= tol)
    cat("   EM iteration:", k, ", Change in 1dfdr fit = ", round(diff, 5), "\n")
    p.0 <- new.p.0; mu.0 <- new.mu.0; sig.0 <- new.sig.0
    f1.tilde <- new.f1.tilde; gam <- new.gam; f <- new.f
  }

  which.z <- gam <= thre.localFDR
  thre <- min(z[which.z])

  if (doplot) {
    hist(z,
         nclass = max(round(length(z)/20), 24),
         probability = TRUE,
         col = "gray", border = "white",
         xlab = "",
         main = "",
         sub = substitute(
           paste(p[0], " = ", p0, ", ",
                 mu[0], " = ", mu0, ", ",
                 sigma[0], " = ", sigma0, ", ",
                 "threshold = ", threshold,
                 sep = ""),
           list(p0 = round(p.0, 2),
                mu0 = round(mu.0, digits = 2),
                sigma0 = round(sig.0, digits = 2),
                threshold = round(thre, digits = 2))))
    rug(z, col = "gray")
    rug(z[which.z], col = 2)
    zs <- sort(z)
    lines(zs, p.0*dnorm(zs, mean = mu.0, sd = sig.0), col = 3, lwd = 2)
    lines(zs, (1-p.0)*f1.tilde[order(z)], col = 2, lwd = 2)
    points(thre, 0, bg = "yellow", col = 2, pch = 25)
  }

  res <- list(p.0 = p.0, mu.0 = mu.0, sig.0 = sig.0, f = f,
              localfdr = gam, iter = k)

  return(res)
}

library(ggpubr)
library(cowplot)


  #library(LogConcDEAD)
  library(logcondens)
  
  z <- as.numeric(z)
  n <- length(z)
  
  ## Initial step
  q0 <- quantile(z, probs = .9)
  p.0 <- mean(z >= q0)
  mu.0 <- mean(z[z >= q0])
  sig.0 <- sd(z[z >= q0])
  
  mu.1 <- mean(z[z < q0])
  sig.1 <- sd(z[z < q0])
  f1.tilde <- dnorm(z, mean = mu.1, sd = sig.1)
  
  f <- gam <- rep(0, n)
  
  ## EM-step
  k <- 0; converged <- 0
  while ( (k < 3) | ((k < max.iter) & (!converged)) ) {
    k <- k + 1
    ## E-step
    tmp <- p.0*dnorm(z, mu.0, sig.0)
    new.f <- tmp + (1 - p.0)*f1.tilde
    new.gam <- tmp/new.f
    
    ## M-step
    w.gam <- new.gam/sum(new.gam, na.rm = TRUE)
    new.mu.0 <- sum(w.gam*z, na.rm = TRUE)
    new.sig.0 <- sqrt(sum(w.gam*(z-new.mu.0)^2, na.rm = TRUE))
    new.p.0 <- mean(new.gam, na.rm = TRUE)
    
    new.f1.tilde <- rep(0, n)
    which.z <- new.gam <= thre.z
    weight <- 1 - new.gam[which.z]
    weight <- weight/sum(weight)
    lcd <- mlelcd(x = z[which.z], w = weight)
    new.f1.tilde[which.z] <- exp(lcd$logMLE)[rank(z[which.z])]
    which.gam <- (new.gam <= Uthre.gam)*(new.gam >= Lthre.gam)
    diff <- max(abs(gam - new.gam)[which.gam])
    converged <- (diff <= tol)
    cat("   EM iteration:", k, ", Change in 1dfdr fit = ", round(diff, 5), "\n")
    p.0 <- new.p.0; mu.0 <- new.mu.0; sig.0 <- new.sig.0
    f1.tilde <- new.f1.tilde; gam <- new.gam; f <- new.f
  }
  
  which.z <- gam <= thre.localFDR
  thre <- max(z[which.z])
  
  if (doplot) {
    hist(z,
         nclass = 100,
         probability = TRUE,
         col = "gray", border = "white",
         xlab = "",
         main = "",
         sub = substitute(
           paste(p[0], " = ", p0, ", ",
                 mu[0], " = ", mu0, ", ",
                 sigma[0], " = ", sigma0, ", ",
                 "threshold = ", threshold,
                 sep = ""),
           list(p0 = round(p.0, 2),
                mu0 = round(mu.0, digits = 2),
                sigma0 = round(sig.0, digits = 2),
                threshold = round(thre, digits = 2))))
    rug(z, col = "#00AFBB")
    rug(z[which.z], col = "#E7B800")
    zs <- sort(z)
    
    lines(zs, p.0*dnorm(zs, mean = mu.0, sd = sig.0), col = "#00AFBB", lwd = 2)
    lines(zs, (1-p.0)*f1.tilde[order(z)], col = "#E7B800", lwd = 2)
    points(thre, 0, bg = "yellow", col = "#E7B800", pch = 25)
    abline(v=mu.0,col="#00AFBB",lty='dashed')
  }
  
  res <- list(p.0 = p.0, mu.0 = mu.0, sig.0 = sig.0, f = f,
              localfdr = gam, iter = k)


sp.mix.1D(z_vel, tol = 5.0e-6, max.iter = 30, doplot = TRUE, thre.localFDR = 0.2, thre.z = 0.95, Uthre.gam = 0.9, Lthre.gam = 0.01)

sp.mix.new(z_vel, tol = 5.0e-6, max.iter = 30, doplot = TRUE, thre.localFDR = 0.2, thre.z = 0.95, Uthre.gam = 0.9, Lthre.gam = 0.01)

z=z_vel
tol = 5.0e-6
max.iter = 30
doplot = TRUE
thre.localFDR = 0.2
thre.z = 0.95
Uthre.gam = 0.9
Lthre.gam = 0.01
p.0 = left$p0
gam=left$localFDR
mu.0 = left$mu0
sig.0 = left$sig0
f1.tilde = left$f1
library(ggplot2)

df=data.frame(z=z,)

ggplot(df) + 
  geom_histogram(breaks=10,aes(x=vector,y=..density..), position="identity") + 
  geom_density(aes(x=vector,y=..density..))

ggplot(df,aes(x=z)) + 
  geom_histogram(aes(y = ..density..),colour = 1, fill = "white",bins=100) +
  geom_line(aes(sort(z), p.0*dnorm(zs, mean = mu.0, sd = sig.0)),color = "#E7B800",lwd=1.2) +
  geom_line(aes(sort(z), (1-p.0)*f1.tilde[order(z)]),color = "#00AFBB",lwd=1.2) +
  geom_vline(aes(xintercept=mu.0), color="#E7B800",linetype="dashed") +
  labs(x="z-value", y = "density")+
  ggtitle(sub)+
  geom_rug(aes(z), color = "#00AFBB") +
  geom_rug(aes(z[which.z]), color = "#E7B800") +
  geom_rug(alpha=0.7, position='jitter')+
  theme_minimal() 
  
sub=substitute(
  paste(p[0], " = ", p0, ", ",
        mu[0], " = ", mu0, ", ",
        sigma[0], " = ", sigma0, ", ",
        "threshold = ", threshold,
        sep = ""),
  list(p0 = round(p.0, 2),
       mu0 = round(mu.0, digits = 2),
       sigma0 = round(sig.0, digits = 2),
       threshold = round(thre, digits = 2)))


