devtools::install_github("JaneeChoi/SpMix",ref = "fmlogcondens_merge")
library(SpMix)
set.seed(210828)
source("/Users/user/Documents/GitHub/Semiparam_Mixture/test_case/fmlogcondens_test/fm_functions.R")
setwd("/Users/user/Documents/GitHub/Semiparam_Mixture/test_case/fmlogcondens_test/Results_faster")

for (N in c(100,500)){
  Res.1<-SimMultNormal(M = 100, n = N, p0 = 0.95)
  Res.2<-SimMultNormal(M = 100, n = N, p0 = 0.90)
  Res.3<-SimMultNormal(M = 100, n = N, p0 = 0.80)
  result_df_1<-rbind(Res.1,Res.2,Res.3)
  write.csv(result_df_1,paste0("2D_Normal_",N,".csv"))
}



for (N in c(100,500)){
  Res.1<-SimMultNormal3d(M = 100, n = N, p0 = 0.95)
  Res.2<-SimMultNormal3d(M = 100, n = N, p0 = 0.90)
  Res.3<-SimMultNormal3d(M = 100, n = N, p0 = 0.80)
  result_df_3<-rbind(Res.1,Res.2,Res.3)  
  write.csv(result_df_3,paste0("3D_Normal_",N,".csv"))
}

for (N in c(100,500)){
  Res.1<-SimMultGamma3d(M = 100, n = N, p0 = 0.95)
  Res.2<-SimMultGamma3d(M = 100, n = N, p0 = 0.90)
  Res.3<-SimMultGamma3d(M = 100, n = N, p0 = 0.80)
  result_df_4<-rbind(Res.1,Res.2,Res.3)  
  write.csv(result_df_4,paste0("3D_Normal_",N,".csv"))
}


for (N in c(100,500)){
  Res.1<-SimMultGamma(M = 100, n = N, p0 = 0.95)
  Res.2<-SimMultGamma(M = 100, n = N, p0 = 0.90)
  Res.3<-SimMultGamma(M = 100, n = N, p0 = 0.80)
  result_df_2<-rbind(Res.1,Res.2,Res.3)
  write.csv(result_df_2,paste0("2D_Gamma_",N,".csv"))
}
