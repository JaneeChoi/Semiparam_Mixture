
# installation
library(devtools)
install_github("JungiinChoi/multiLocalFDR")

library(multiLocalFDR)
?multiLocalFDR

# Case Study: Pathways
## Data
data("pathways", package = "multiLocalFDR")
head(pathways)

## -Inf fix
pathways[238,4] = 9.99999e-01

## SpMix for each pathway
params_L <- SpMix(pathways[,2], p_value = TRUE, tol = 1e-10)
params_O <- SpMix(pathways[,3], p_value = TRUE, tol = 1e-10)
params_S <- SpMix(pathways[,4], p_value = TRUE, tol = 1e-10)

## Plot results for each pathway

plotFDR(pathways[,2], params_L$p0, params_L$mu0, params_L$sig0, params_L$f1, 
        params_L$localFDR, p_value = TRUE)
plotFDR(pathways[,3], params_O$p0, params_O$mu0, params_O$sig0, params_O$f1, 
        params_O$localFDR, p_value = TRUE)
plotFDR(pathways[,4], params_S$p0, params_S$mu0, params_S$sig0, params_S$f1, 
        params_S$localFDR, p_value = TRUE)

## SpMix for 3-dimensional data
params <- SpMix(pathways[,2:4], p_value = TRUE)

## Plot results for 3-dimensional data
## To-do
plotFDR(z_galaxy, params$p0, params$mu0, params$sig0, params$f1, 
        params$localFDR, alternative = "less", thre_localFDR = 0.2, testing = FALSE)

# Case Study: Galaxy
## Data
data("galaxy", package = "multiLocalFDR")
head(galaxy)

## transform to z-values
z_galaxy <- galaxy$velocity
head(z_galaxy)

## SpMix for galaxy data
params <- SpMix(z_galaxy, alternative = "less")

## Plot results
plotFDR(z_galaxy, params$p0, params$mu0, params$sig0, params$f1, 
        params$localFDR, alternative = "less", thre_localFDR = 0.2, testing = FALSE)
