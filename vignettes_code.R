# installation
library(devtools)
install_github("JungiinChoi/multiLocalFDR")

library(multiLocalFDR)

# Simulations

## 1d Normal + Normal p0=0.8, N(0,1), N(3.5,0.5^2)
n=5000
p0=0.8

n0 <- rbinom(1, n, p0)
n1 <- n - n0
z0 <- rnorm(n0)
z1 <- rnorm(n1, mean = 3.5, sd = 0.5)
z_NN1d <- c(z0, z1)
NN1d <- SPMix(z_NN1d, thre_z = 0.99, Uthre_gam = 0.99)

plotSPMix(z_NN1d, NN1d$p0, NN1d$mu0, NN1d$sig0, NN1d$f, NN1d$f1, NN1d$localFDR)
plotSPMix(z_NN1d, NN1d$p0, NN1d$mu0, NN1d$sig0, NN1d$f, NN1d$f1, NN1d$localFDR,
          testing = FALSE)

## 1d Normal + Gamma p0=0.8, N(0,1), gamma(12,0.25)

z1_gamma <- rgamma(n1, shape = 12, rate = (1/0.25))
z_NG1d <- c(z0, z1_gamma)
NG1d <- SPMix(z_NG1d, thre_z = 0.99, Uthre_gam = 0.99)
plotSPMix(z_NG1d, NG1d$p0, NG1d$mu0, NG1d$sig0, NG1d$f, NG1d$f1, NG1d$localFDR)
plotSPMix(z_NG1d, NG1d$p0, NG1d$mu0, NG1d$sig0, NG1d$f, NG1d$f1, NG1d$localFDR,
          testing = FALSE)

## 2d Normal + Normal

library(copula)
library(mclust)

n0 <- rbinom(1, n, p0)
n1 <- n - n0
sig0 <- matrix(c(1, 0.25, 0.25, 1), ncol = 2, byrow = T)
V <- rCopula(n1, ellipCopula(family = "normal", param = 0.5, dim = 2))
z0 <- rmvnorm(n0, sigma = sig0)
z1 <- qnorm(V, mean = 3.5, sd = .5)
z_NN2d <- data.frame(rbind(z0, z1))

ggplot(z_NN, aes(x = z_NN[,1], y = z_NN[,2])) +
  geom_point(aes(color = factor(c(rep(1,n0), rep(2,n1)), labels = c("Normal(Null)", "Normal(Alternative)")))) +
  scale_color_manual(values = c("#999999", "#E69F00"), name="")  +
  labs(title = "Random Generated Multi-Normal Mixture", x="x", y = "y") +
  theme_classic()

NN <- SPMix(z_NN, tol = 1e-10)

plotSPMix(z_NN, NN$p0, NN$mu0, NN$sig0, NN$f, NN$f1, NN$localFDR, 
          xlab = "x", ylab = "y")

## 2d Normal + Gamma

z1_gamma <- qgamma(V, shape = 12, rate = 4)
z_NG <- data.frame(rbind(z0, z1_gamma))

ggplot(z_NG, aes(x = z_NG[,1], y = z_NG[,2])) +
  geom_point(aes(color = factor(c(rep(1,n0), rep(2,n1)), labels = c("Normal", "Gamma")))) +
  scale_color_manual(values = c("#999999", "#E69F00"), name="")  +
  labs(title = "Random Generated Normal/Gamma Mixture", x="x", y = "y") +
  theme_classic()

NG <- SPMix(z_NG, tol = 1e-10)
plotSPMix(z_NG, NN$p0, NN$mu0, NN$sig0, NN$f, NN$f1, NN$localFDR, 
          xlab = "x", ylab = "y")

## 3d Normal + Normal

library(scatterplot3d)

Sigma0 <- matrix(c(1, 0.25, 0.25, 0.25, 1, 0.25, 0.25,0.25,1), ncol = 3)
n0 <- rbinom(1, n, p0)
n1 <- n - n0
V <- rCopula(n1, ellipCopula(family = "normal", param = 0.5, dim = 3))
z0_NN3d <- rmvnorm(n0, sigma = Sigma0)
z1_NN3d <- qnorm(V, mean = 3.5, sd = .5)
z_NN3d <- rbind(z0_NN3d, z1_NN3d)

NN3d <- SPMix(z_NN3d, 1e-10)

colors <- c("#999999", "#E69F00")
colors <- colors[as.numeric(NN3d$localFDR <= 0.1)+1]
scatterplot_NN<- scatterplot3d(z_NN3d, pch = 16, color=colors,
                            xlab = "x", ylab = "y", zlab = "z")
legend_pathway <- factor((NN3d$localFDR <= 0.1), levels = c(FALSE, TRUE), 
                         labels = c("Nonsignificant", "Significant"))
legend(scatterplot_NN$xyz.convert(6, -5, 0), legend = levels(legend_pathway),
       col = c("#999999", "#E69F00"), pch = 16)

## 3d Normal + Gamma

V_NG <- rCopula(n1, frankCopula(param = param, dim = 3))
z0_NG3d <- rmvnorm(n0, sigma = Sigma0)
z1_NG3d <- qgamma(V_NG, shape = 12, rate = 4)
z_NG3d <- rbind(z0_NG3d, z1_NG3d)

NG3d <- SPMix(z_NG3d, 1e-10)

colors <- c("#999999", "#E69F00")
colors <- colors[as.numeric(NG3d$localFDR <= 0.1)+1]
scatterplot_NG<- scatterplot3d(z_NG3d, pch = 16, color=colors,
                            xlab = "x", ylab = "y", zlab = "z")
legend_pathway <- factor((NG3d$localFDR <= 0.2), levels = c(FALSE, TRUE), 
                         labels = c("Nonsignificant", "Significant"))
legend(scatterplot_NG$xyz.convert(6, -5, 0), legend = levels(legend_pathway),
       col = c("#999999", "#E69F00"), pch = 16)

# Case Study: Pathways
## Data
data("pathways", package = "multiLocalFDR")
?pathways
head(pathways)

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

## SpMix for 3-dimensional data
params_LOS <- SPMix(pathways[,2:4], p_value = TRUE, tol = 1e-10)

## Pathways which md-fdr <= 0.2
head(pathways[params_LOS$localFDR <= 0.1,]$Gene.Set)
length(pathways[params_LOS$localFDR <= 0.1,]$Gene.Set)

# md-fdr dataframe
pathways_df <- pathways
pathways_df$localFDR <- params_LOS$localFDR
pathways_df$sgnf <- (params_LOS$localFDR <= 0.01)
head(pathways_df)

save_df <-pathways_df[pathways_df$sgnf, ]
save_df <- save_df[order(save_df$localFDR),]
save_df


## 3d scatter plot
install.packages("scatterplot3d") # Install
library(scatterplot3d)

colors <- c("#999999", "#E69F00")
colors <- colors[as.numeric(pathways_df$sgnf)+1]
scatterplot<- scatterplot3d(qnorm(1-as.matrix(pathways_df[,2:4])), pch = 16, color=colors,
              xlab = "Peripheral Leukocytes",
              ylab = "Orbital Tissue",
              zlab = "Sinus Brushings")

legend_pathway <- factor(pathways_df$sgnf, levels = c(FALSE, TRUE), labels = c("Nonsignificant", "Significant"))
legend(scatterplot$xyz.convert(6, -5, 0), legend = levels(legend_pathway),
       col = c("#999999", "#E69F00"), pch = 16)

## Animated 3d Plot
### install.packages(c("rgl","magick"))
library( rgl )
plot3d(qnorm(1-as.matrix(pathways_df[,2:4])), col = colors, type = "s", radius = .2)
play3d( spin3d( axis = c(0, 0, 1), rpm = 20), duration = 10 )
movie3d(
  movie="3dAnimatedScatterplot", spin3d( axis = c(0, 0, 1), rpm = 7),
  duration = 10, dir = "/Users/user/Desktop", type = "gif", clean = TRUE)


# Case Study: Galaxy
## Data
data("galaxy", package = "multiLocalFDR")
?galaxy
head(galaxy)

## transform to z-values
z_galaxy <- galaxy$velocity
head(z_galaxy)

## SpMix for galaxy data
params <- SPMix(z_galaxy, alternative = "less")
params

## Plot results
plotSPMix(z_galaxy, params$p0, params$mu0, params$sig0, params$f, params$f1, 
          params$localFDR, alternative = "greater", testing = FALSE,
                      xlab = "velocity(km/s)")
