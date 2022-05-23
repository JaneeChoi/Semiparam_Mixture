# Load Library
library(LogConcDEAD)
library(mclust)
library(logcondens)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(grid)

# Load Data
rv<-read.table("rvdata_for_stat.txt", skip = 1)
z_vel = scale(rv$V2)[,1]

sp.mix.1D <- function(z, tol = 5.0e-6, max.iter = 30, doplot = TRUE, thre.localFDR = 0.2, thre.z = 0.95, Uthre.gam = 0.9, Lthre.gam = 0.01)
{

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

  res <- list(p.0 = p.0, mu.0 = mu.0, sig.0 = sig.0, f = f,
              localfdr = gam, iter = k)

  return(res)
}


gal_res=sp.mix.1D(z_vel, tol = 5.0e-6, max.iter = 5, doplot = FALSE, thre.localFDR = 0.2, thre.z = 0.95, Uthre.gam = 0.9, Lthre.gam = 0.01)

gal_res = res
z=z_vel
p.0 = gal_res$p0
gam=gal_res$localFDR
mu.0 = gal_res$mu0
sig.0 = gal_res$sig0
f1.tilde = gal_res$f1
which.z = (gam <= 0.2)
thre=min(z[which.z])

n=length(z)
distribution=rep("",n)
for (k in 1:n){
  if (which.z[k]) {
    distribution[k] = "Alternative"
  }
  else{
    distribution[k] = "Null"
  }
}

distribution <- factor(distribution, levels = c("Null","Alternative"))

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

df = data.frame(z=z)
zs <- sort(z)
ggplot(df,aes(x=z)) + 
  geom_histogram(aes(y = ..density..),colour = 1, fill = "white",bins=100) +
  geom_line(aes(sort(z), p.0*dnorm(zs, mean = mu.0, sd = sig.0)),color = "#00BFC4",lwd=1.1) +
  geom_line(aes(sort(z), ((1-p.0)*f1.tilde[order(z)])),color = "#F8766D",lwd=1.1) +
  geom_vline(aes(xintercept=mu.0), color="#00BFC4",linetype="dashed") +
  geom_point(mapping = aes(x = thre, y = 0.01),size = 2,color='yellow',shape=25,fill="yellow") +
  labs(x="z-value", y = "density") +
  ggtitle(sub) +
  theme(plot.title = element_text(margin = margin(b = -10))) +
  geom_rug(aes(z,color = distribution))+
  scale_color_manual(values = c("#00BFC4", "#F8766D"),name="")+
  theme_classic() 

