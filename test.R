library("myPackage")

setwd("Documents/GitHub/Semiparam_Mixture/SPMix-LocalFDR/Carina")
dat <- read.table('carina.dat')
str(dat)

x <- dat[dat$V8 + dat$V9 > 0,]
x <- x[x$V6 < 3,]
vel <- x$V4 # represents the Radial velocity data of stars in the Carina galaxy
vel
vel.fit <- sp.mix.1D(vel)

??myPackage


