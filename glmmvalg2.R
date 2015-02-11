library(GLMMF)

load("borneodata.rda")
load("initmat.rda")
set.seed(12345)
randu1<-rnorm(37)
randu2<-rnorm(37)

fit_glmmf2<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
                  convtol=1e-10,maxiter=25,maxiter2=0,init.factor=initmat_glmmf,
                  data=borneodata,control=list(trace=1,REPORT=1,maxit=2),method="BFGS",
                  estimate=T,trace=1)
