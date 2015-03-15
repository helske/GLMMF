library(GLMMF)
load("spiderdata.rda")
set.seed(1)
randu1<-rnorm(28)
randu2<-c(0,rnorm(27))
fitinit<-glmmf(group="species",response="abundance",distinct=~randu1+randu2,nfactors=0,distribution="poisson",
               data=spiderdata,control=list(trace=1,REPORT=1),trace=2,estimate=F,method="BFGS")
fitinit$log #-3436.87
initmat<-t(matrix(fitinit$coefs[,28],nrow=3)[2:3,])
fit<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
           data=spiderdata,control=list(trace=1,REPORT=1,maxit=10000),
           method="BFGS",init.factor=initmat,estimate=T,trace=1,convtol=1e-7)

fit$log #857.248226 
factors<-fit$coefs[1:2,]
plot(x=factors[1,],y=factors[2,],pch=NA)
text(x=factors[1,],y=factors[2,], labels=as.character(1:28))


a<-proc.time()
fit1<-optim(par=model$Z[,13:14,1][lower.tri(matrix(0,12,2),TRUE)],
            fn=likfn, model=model,method='BFGS')
proc.time()-a

a<-proc.time()
fit<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
           data=spiderdata,
           method="BFGS",init.factor=initmat,estimate=T,trace=0,convtol=1e-7)
proc.time()-a

