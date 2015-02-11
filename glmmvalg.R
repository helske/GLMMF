load("borneodata.rda")
library(GLMMF)

set.seed(12345)
randu1<-rnorm(37)
randu2<-rnorm(37)
fitinit<-glmmf(group="species",response="abundance",distinct=~Logging+Slope+randu1+randu2,
               nfactors=0,distribution="poisson",maxiter=100,
               data=borneodata,control=list(trace=1,REPORT=1),method="BFGS",trace=2)

fitinit$log #-1433.71
initmat<-t(matrix(fitinit$coefs[,37],nrow=7)[6:7,])

fit<-glmmf(group="species",response="abundance",distinct=~Logging+Slope,nfactors=2,distribution="poisson",
           convtol=1e-8,maxiter=10,maxiter2=0,init.theta=fitinit$linear.predictor,
           data=borneodata,control=list(trace=1,REPORT=1,maxit=5),method="BFGS",init.factor=initmat,estimate=T,trace=0)
fit$log
