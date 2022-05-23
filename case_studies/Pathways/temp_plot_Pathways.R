# Load Library
library(LogConcDEAD)
library(mclust)
library(logcondens)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(grid)


# Load Data
load("Study1.RData")
z1 <- read.csv("L.csv")
z2 <- read.csv("O.csv")
z3 <- read.csv("S.csv")
localfdr_df <- read.csv("localFDR.csv")

gs <- unique(c(as.character(z1[,1]),
               as.character(z2[,1]),
               as.character(z3[,1])))

id1 <- match(gs, z1[,1])
id2 <- match(gs, z2[,1])
id3 <- match(gs, z3[,1])

id <- cbind(id1, id2, id3)

tmp <- apply(id, 1, is.na)
tmp <- apply(tmp, 2, sum)
idA <- id[!is.na(id[,1]) & !is.na(id[,2])&!is.na(id[,3]),]

z <- cbind(z1[idA[,1],], z2[idA[,2],], z3[idA[,3],])

z0 <- z[,c(5, 12, 19)]
z0 <- -qnorm(as.matrix(z0))
z0[,3] <- pmax(z0[,3], -3)
res_list = list(res1,res2,res3)
localfdr_df[,30]


# Plot

par(mfrow=c(3,1))

i=3

for (i in 1:3){
  which.z=localfdr_df[,(29+i)]
  z = z0[,i]
  thre <- min(z[which.z])
  res = res_list[[i]]
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
    geom_line(aes(sort(z), dnorm(zs, mean = mu.0, sd = sig.0)/p.0),color = "#00BFC4",lwd=1.1) +
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



