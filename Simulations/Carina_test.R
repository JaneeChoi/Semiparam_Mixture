devtools::install_github("JaneeChoi/SpMix")
library(SpMix)

?sp.mix.1D

dat <- read.table('/Users/user/Downloads/SPMix-LocalFDR-master/Carina/carina.dat')
str(dat)

x <- dat[dat$V8 + dat$V9 > 0,]
x <- x[x$V6 < 3,]
vel <- x$V4 # represents the Radial velocity data of stars in the Carina galaxy
vel
Carinadata<-data.frame(vel=vel)
head(Carinadata)
vel

?sp.mix.1D

vel.fit <- sp.mix.1D(vel,doplot=FALSE)
