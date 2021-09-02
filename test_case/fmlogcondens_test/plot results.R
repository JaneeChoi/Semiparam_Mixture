
# LogConcDEAD

# normal

setwd("/Users/choiiiiii/Documents/GitHub/Semiparam_Mixture/test_case/fmlogcondens_test/Multi_Gamma_LocConcDEAD")


res100<-unname(as.matrix(read.csv("result_gamma_100.csv")[,2:4]))
res500<-unname(as.matrix(read.csv("result_gamma_500.csv")[,2:4]))
res1000<-unname(as.matrix(read.csv("result_gamma_1000.csv")[,2:4]))
res5000<-unname(as.matrix(read.csv("result_gamma_5000.csv")[,2:4]))
res10000<-unname(as.matrix(read.csv("result_gamma_10000.csv")[,2:4]))
#res15000<-unname(as.matrix(read.csv("result_15000.csv")[,2:4]))
#res20000<-unname(as.matrix(read.csv("result_20000.csv")[,2:4]))
#res25000<-unname(as.matrix(read.csv("result_25000.csv")[,2:4]))


normal.log<-data.frame()
normal.log<-rbind(normal.log,data.frame(N=rep(100,3),p0=c(0.95,0.9,0.8),p0_hat=res100[2,],time=res100[1,],sens=res100[3,],spec=res100[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(500,3),p0=c(0.95,0.9,0.8),p0_hat=res500[2,],time=res500[1,],sens=res500[3,],spec=res500[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(1000,3),p0=c(0.95,0.9,0.8),p0_hat=res1000[2,],time=res1000[1,],sens=res1000[3,],spec=res1000[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(5000,3),p0=c(0.95,0.9,0.8),p0_hat=res5000[2,],time=res5000[1,],sens=res5000[3,],spec=res5000[4,]))
normal.log<-rbind(normal.log,data.frame(N=rep(10000,3),p0=c(0.95,0.9,0.8),p0_hat=res10000[2,],time=res10000[1,],sens=res10000[3,],spec=res10000[4,]))
#normal.log<-rbind(normal.log,data.frame(N=rep(15000,3),p0=c(0.95,0.9,0.8),p0_hat=res15000[2,],time=res15000[1,],sens=res15000[3,],spec=res15000[4,]))
#normal.log<-rbind(normal.log,data.frame(N=rep(20000,3),p0=c(0.95,0.9,0.8),p0_hat=res20000[2,],time=res20000[1,],sens=res20000[3,],spec=res20000[4,]))
#normal.log<-rbind(normal.log,data.frame(N=rep(25000,3),p0=c(0.95,0.9,0.8),p0_hat=res25000[2,],time=res25000[1,],sens=res25000[3,],spec=res25000[4,]))

write.csv(normal.log,"merged.csv")

library(ggplot2)
library(dplyr)

p1<-ggplot(data=normal.log,aes(x=N,y=time,group=p0,colour=p0))+geom_line()
ggsave("time.png", plot=p1, height=6, width=8, dpi=600)

p2<-ggplot(data=normal.log,aes(x=N,y=sens,group=p0,colour=p0))+geom_line()+coord_cartesian(ylim = c(0.5,1)) +ylab("sensitivity")
p2
ggsave("sens.png", plot=p2, height=6, width=8, dpi=600)

p3<-ggplot(data=normal.log,aes(x=N,y=spec,group=p0,colour=p0))+geom_line()+coord_cartesian(ylim = c(0.5,1)) +ylab("specificity")
p3
ggsave("spec.png", plot=p3, height=6, width=8, dpi=600)

p4<-ggplot(data=normal.log,aes(x=N,y=p0_hat,group=p0,colour=p0))+coord_cartesian(xlim=c(50,10000))+geom_line()+ylab("specificity")
p4
ggsave("p0_hat.png", plot=p4, height=6, width=8, dpi=600)


# gamma

setwd("/Users/choiiiiii/Documents/GitHub/Semiparam_Mixture/test_case/fmlogcondens_test/Multi_Gamma_LocConcDEAD")

result_gamma_100<-unname(as.matrix(read.csv("result_gamma_100.csv")[,2:4]))
result_gamma_500<-unname(as.matrix(read.csv("result_gamma_500.csv")[,2:4]))
result_gamma_1000<-unname(as.matrix(read.csv("result_gamma_1000.csv")[,2:4]))
result_gamma_5000<-unname(as.matrix(read.csv("result_gamma_5000.csv")[,2:4]))
result_gamma_10000<-unname(as.matrix(read.csv("result_gamma_10000.csv")[,2:4]))


gamma.log<-data.frame()
gamma.log<-rbind(gamma.log,data.frame(N=rep(100,3),p0=c(0.95,0.9,0.8),p0_hat=result_gamma_100[2,],time=result_gamma_100[1,],sens=result_gamma_100[3,],spec=result_gamma_100[4,]))
gamma.log<-rbind(gamma.log,data.frame(N=rep(500,3),p0=c(0.95,0.9,0.8),p0_hat=result_gamma_500[2,],time=result_gamma_500[1,],sens=result_gamma_500[3,],spec=result_gamma_500[4,]))
gamma.log<-rbind(gamma.log,data.frame(N=rep(1000,3),p0=c(0.95,0.9,0.8),p0_hat=result_gamma_1000[2,],time=result_gamma_1000[1,],sens=result_gamma_1000[3,],spec=result_gamma_1000[4,]))
gamma.log<-rbind(gamma.log,data.frame(N=rep(5000,3),p0=c(0.95,0.9,0.8),p0_hat=result_gamma_5000[2,],time=result_gamma_5000[1,],sens=result_gamma_5000[3,],spec=result_gamma_5000[4,]))
gamma.log<-rbind(gamma.log,data.frame(N=rep(10000,3),p0=c(0.95,0.9,0.8),p0_hat=result_gamma_10000[2,],time=result_gamma_10000[1,],sens=result_gamma_10000[3,],spec=result_gamma_10000[4,]))

write.csv(gamma.log,"merged.csv")

gamma.log
p1<-ggplot(data=gamma.log,aes(x=N,y=time,group=p0,colour=p0))+geom_line()
ggsave("time.png", plot=p1, height=6, width=8, dpi=600)

p2<-ggplot(data=gamma.log,aes(x=N,y=sens,group=p0,colour=p0))+geom_line()+coord_cartesian(ylim = c(0.5,1)) +ylab("sensitivity")
p2
ggsave("sens.png", plot=p1, height=6, width=8, dpi=600)

p3<-ggplot(data=gamma.log,aes(x=N,y=spec,group=p0,colour=p0))+geom_line()+coord_cartesian(ylim = c(0.9,1)) +ylab("specificity")
p3
ggsave("spec.png", plot=p1, height=6, width=8, dpi=600)



