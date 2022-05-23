# Load Library
library(LogConcDEAD)
library(mclust)
library(logcondens)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(plot3D)


# Load Data
load("Study2.RData")
z0=z
res_list = list(Res.1D.1,Res.1D.2)

# Plot
i=2

for (i in 1:3){
  res = res_list[[i]]
  z = z0[,i]
  which.z=res$localfdr <= 0.05
  thre <- min(z[which.z])
  p.0 = res$p.0
  gam=res$localfdr
  mu.0 = res$mu.0
  sig.0 = res$sig.0
  f1.tilde = res$f
  
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
    geom_line(aes(sort(z), dnorm(zs, mean = mu.0, sd = sig.0)),color = "#00BFC4",lwd=1.1) +
    geom_line(aes(sort(z), (f1.tilde[order(z)]-p.0*dnorm(zs, mean = mu.0, sd = sig.0))/0.6),color = "#F8766D",lwd=1.1) +
    geom_vline(aes(xintercept=mu.0), color="#00BFC4",linetype="dashed") +
    geom_point(mapping = aes(x = thre, y = 0.01),size = 2,color='yellow',shape=25,fill="yellow") +
    labs(x="z-value", y = "density")+
    ggtitle(sub) +
    theme(plot.title = element_text(margin = margin(b = -10)))+
    geom_rug(aes(z,color = distribution))+
    scale_color_manual(values = c("#00BFC4", "#F8766D"),name="")+
    theme_classic() 
}

## 2D plot

n=length(z0[,1])
col2d=rep("",n)
for (k in 1:n){
  if (Res.2D$localfdr[k] <= .05) {
    col2d[k] = "#F8766D"
  }
  else{
    col2d[k] = "#00BFC4"
  }
}

n=length(z0[,1])
colrug=rep("",n)
for (k in 1:n){
  if (Res.2D$localfdr[k] <= .05) {
    colrug[k] = "Alternative"
  }
  else{
    colrug[k] = "Null"
  }
}
colrug <- factor(colrug, levels = c("Null","Alternative"))


ggplot() +
  geom_point(mapping = aes(x = z0[,1], y = z0[,2]), color = col2d) +
  labs(x="z1", y = "z2") +
  geom_rug(mapping = aes(x = z0[,1], y = z0[,2], color = colrug)) +
  scale_color_manual(values = c("#00BFC4", "#F8766D"), name="") +
  theme_classic() 

## 3D plot

scatter3D_fancy <- function(x, y, z,..., colvar = z)
{
  panelfirst <- function(pmat) {
    XY <- trans3D(x, y, z = rep(min(z), length(z)), pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = colvar, pch = ".", 
              cex = 2, add = TRUE, colkey = FALSE)
    
    XY <- trans3D(x = rep(min(x), length(x)), y, z, pmat = pmat)
    scatter2D(XY$x, XY$y, colvar = colvar, pch = ".", 
              cex = 2, add = TRUE, colkey = FALSE)
  }
  scatter3D(x, y, z, ..., colvar = colvar, panel.first=panelfirst,
            colkey = list(length = 0.5, width = 0.5, cex.clab = 0.75)) 
}

n=length(z0[,1])
col2d=rep("",n)
for (k in 1:n){
  if (Res.2D$localfdr[k] <= .05) {
    col2d[k] = "#F8766D"
  }
  else{
    col2d[k] = "#00BFC4"
  }
}

scatter3D_fancy(z0[,1], z0[,2], Res.2D$localfdr, pch = ".", bty="g",
                main = "Iris data",  colvar = as.numeric(Res.2D$localfdr <= 0.05))


as.numeric(Res.2D$localfdr <= 0.05)

