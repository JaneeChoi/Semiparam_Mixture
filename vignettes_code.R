# Installation
library(devtools)
install_github("JungiinChoi/multiLocalFDR")

## required libraries
library(copula)
library(mvtnorm)
library(ggplot2)
library(scatterplot3d)

library(multiLocalFDR)

# Simulations

## Univaraite

### Normal + Normal

# set seed
set.seed(1)

n <- 10000
p0 <- 0.8
n0 <- rbinom(1, n, p0)
n1 <- n - n0
z0 <- rnorm(n0)
z1 <- rnorm(n1, mean = 3.5, sd = 0.5)
z_NN1d <- c(z0, z1)

NN1d <- SPMix(z_NN1d)
str(NN1d)

# Hypothesis testing 
plotSPMix(z_NN1d, NN1d$p0, NN1d$mu0, NN1d$sig0, NN1d$f, NN1d$f1, NN1d$localFDR)

# Density Estimation
plotSPMix(z_NN1d, NN1d$p0, NN1d$mu0, NN1d$sig0, NN1d$f, NN1d$f1, NN1d$localFDR,
          testing = FALSE)

### Normal + Gamma

z1_gamma <- rgamma(n1, shape = 12, rate = 4)
z_NG1d <- c(z0, z1_gamma)

NG1d <- SPMix(z_NG1d)

# Hypothesis testing 
plotSPMix(z_NG1d, NG1d$p0, NG1d$mu0, NG1d$sig0, NG1d$f, NG1d$f1, NG1d$localFDR)

# Density Estimation
plotSPMix(z_NG1d, NG1d$p0, NG1d$mu0, NG1d$sig0, NG1d$f, NG1d$f1, NG1d$localFDR,
          testing = FALSE)

## 2-dimensional data

### Normal + Gaussian Copula
n <- 3000
n0 <- rbinom(1, n, p0)
n1 <- n - n0
sig0 <- matrix(c(1, 0.25, 0.25, 1), ncol = 2, byrow = T)
V <- rCopula(n1, ellipCopula(family = "normal", param = 0.5, dim = 2))
z0 <- rmvnorm(n0, sigma = sig0)
z1 <- qnorm(V, mean = 3.5, sd = .5)
z_NN2d <- data.frame(rbind(z0, z1))

ggplot(z_NN2d, aes(x = z_NN2d[,1], y = z_NN2d[,2])) +
  geom_point(aes(color = factor(c(rep(1,n0), rep(2,n1)), labels = c("Normal(Null)", "Gaussian Copula(Alternative)")))) +
  scale_color_manual(values = c("#999999", "#E69F00"), name="")  +
  labs(title = "Random Generated Mixture density data", x="x", y = "y") +
  theme_classic()

NN2d <- SPMix(z_NN2d)

plotSPMix(z_NN2d, NN2d$p0, NN2d$mu0, NN2d$sig0, NN2d$f, NN2d$f1, NN2d$localFDR)

### Normal + Frank Copula

z1_gamma <- qgamma(V, shape = 12, rate = 4)
z_NG2d <- data.frame(rbind(z0, z1_gamma))

ggplot(z_NG2d, aes(x = z_NG2d[,1], y = z_NG2d[,2])) +
  geom_point(aes(color = factor(c(rep(1,n0), rep(2,n1)), labels = c("Normal", "Gamma")))) +
  scale_color_manual(values = c("#999999", "#E69F00"), name="")  +
  labs(title = "Random Generated Normal/Gamma Mixture", x="x", y = "y") +
  theme_classic()

NG2d <- SPMix(z_NG2d)

plotSPMix(z_NG2d, NG2d$p0, NG2d$mu0, NG2d$sig0, NG2d$f, NG2d$f1, NG2d$localFDR)

## 3-dimensional data

### Normal + Gaussian Copula

Sigma0 <- matrix(c(1, 0.25, 0.25, 0.25, 1, 0.25, 0.25,0.25,1), ncol = 3)
n0 <- rbinom(1, n, p0)
n1 <- n - n0
V <- rCopula(n1, ellipCopula(family = "normal", param = 0.5, dim = 3))
z0_NN3d <- rmvnorm(n0, sigma = Sigma0)
z1_NN3d <- qnorm(V, mean = 3.5, sd = .5)
z_NN3d <- rbind(z0_NN3d, z1_NN3d)

NN3d <- SPMix(z_NN3d)

plotSPMix(z_NN3d, NN3d$p0, NN3d$mu0, NN3d$sig0, NN3d$f, NN3d$f1, NN3d$localFDR, thre_localFDR = 0.1)

### Normal + Frank Copula

V_NG <- rCopula(n1, frankCopula(param = 0.5, dim = 3))
z0_NG3d <- rmvnorm(n0, sigma = Sigma0)
z1_NG3d <- qgamma(V_NG, shape = 12, rate = 4)
z_NG3d <- rbind(z0_NG3d, z1_NG3d)

NG3d <- SPMix(z_NG3d)

plotSPMix(z_NG3d, NG3d$p0, NG3d$mu0, NG3d$sig0, NG3d$f, NG3d$f1, NG3d$localFDR, thre_localFDR = 0.1)

# Applications

## Galaxy Velocity (Univariate)

data("galaxy", package = "multiLocalFDR")
head(galaxy)

z_galaxy <- galaxy$velocity
params <- SPMix(z_galaxy, alternative = "less")

plotSPMix(z_galaxy, params$p0, params$mu0, params$sig0, params$f, params$f1, 
          params$localFDR, alternative = "greater", testing = FALSE,
          xlab = "velocity(km/s)")

## Pathways (Multivariate)

data("pathways", package = "multiLocalFDR")
head(pathways)

### Univariate analysis

## SpMix for each pathway
params_L <- SPMix(pathways[,2], p_value = TRUE, tol = 1e-10)
params_O <- SPMix(pathways[,3], p_value = TRUE, tol = 1e-10)
params_S <- SPMix(pathways[,4], p_value = TRUE, tol = 1e-10)

## Plot results for each pathway
plotSPMix(pathways[,2], params_L$p0, params_L$mu0, params_L$sig0, params_L$f, 
          params_L$f1, params_L$localFDR, p_value = TRUE, xlab = "Peripheral Leukocytes")
plotSPMix(pathways[,3], params_O$p0, params_O$mu0, params_O$sig0, params_O$f, 
          params_O$f1, params_O$localFDR, p_value = TRUE, xlab = "Orbital Tissue")
plotSPMix(pathways[,4], params_S$p0, params_S$mu0, params_S$sig0, params_S$f, 
          params_S$f1, params_S$localFDR, p_value = TRUE, xlab = "Sinus Brushings")

### Multivariate Analysis

## SpMix for 3-dimensional data
params_LOS <- SPMix(pathways[,2:4], p_value = TRUE)

## Pathways which md-fdr <= 0.01
head(pathways[params_LOS$localFDR <= 0.01,]$Gene.Set)
length(pathways[params_LOS$localFDR <= 0.01,]$Gene.Set)

### 3d scatter plot 


plotSPMix(pathways[,2:4], params_LOS$p0, params_LOS$mu0, params_LOS$sig0, params_LOS$f, params_LOS$f1,
          params_LOS$localFDR, p_value = TRUE, thre_localFDR = 0.01, xlab = "Peripheral Leukocytes",
          ylab = "Orbital Tissue", zlab = "Sinus Brushings", coord_legend = c(6, -5, 0))
