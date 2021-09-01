setwd("/Users/user/Documents/GitHub/Semiparam_Mixture/test_case/fmlogcondens_test/Multi_Normal_LocConcDEAD")


res100<-unname(as.matrix(read.csv("result_100.csv")[,2:4]))
res500<-unname(as.matrix(read.csv("result_500.csv")[,2:4]))
res1000<-unname(as.matrix(read.csv("result_1000.csv")[,2:4]))
res5000<-unname(as.matrix(read.csv("result_5000.csv")[,2:4]))
res10000<-unname(as.matrix(read.csv("result_10000.csv")[,2:4]))
res15000<-unname(as.matrix(read.csv("result_15000.csv")[,2:4]))
res20000<-unname(as.matrix(read.csv("result_20000.csv")[,2:4]))
res25000<-unname(as.matrix(read.csv("result_25000.csv")[,2:4]))

normal.log<-data.frame()
normal.log<-rbind(normal.log,data.frame(N=rep(100,3),p0=c(0.95,0.9,0.8),p0_hat=res100[2,],time=res100[1,],sens=res100[3,],spec=res100[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(500,3),p0=c(0.95,0.9,0.8),p0_hat=res500[2,],time=res500[1,],sens=res500[3,],spec=res500[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(1000,3),p0=c(0.95,0.9,0.8),p0_hat=res1000[2,],time=res1000[1,],sens=res1000[3,],spec=res1000[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(5000,3),p0=c(0.95,0.9,0.8),p0_hat=res5000[2,],time=res5000[1,],sens=res5000[3,],spec=res5000[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(10000,3),p0=c(0.95,0.9,0.8),p0_hat=res10000[2,],time=res10000[1,],sens=res10000[3,],spec=res10000[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(15000,3),p0=c(0.95,0.9,0.8),p0_hat=res15000[2,],time=res15000[1,],sens=res15000[3,],spec=res15000[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(20000,3),p0=c(0.95,0.9,0.8),p0_hat=res20000[2,],time=res20000[1,],sens=res20000[3,],spec=res20000[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(25000,3),p0=c(0.95,0.9,0.8),p0_hat=res25000[2,],time=res25000[1,],sens=res25000[3,],spec=res25000[4,]))

normal.log
