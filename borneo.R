library(GLMMF)
dat = read.csv("born.csv")
Y = dat[,-(1:4)]
X = dat[,2:4] #remove id column
#remove rare species:
Y = as.matrix(Y[,apply(Y>0,2,sum)>30])

Xdat<-do.call(rbind, replicate(ncol(Y), X, simplify=FALSE)) # 5069  3
species<-rep(colnames(Y),each=nrow(Y))
borneodata<-data.frame(c(Y),Xdat,species)
names(borneodata)<-c("abundance",colnames(Xdat),"species")
save(borneodata,file="borneodata.rda")

load("borneodata.rda")


set.seed(12345)
randu1<-rnorm(37)
randu2<-rnorm(37)
fitinit<-glmmf(group="species",response="abundance",distinct=~randu1+randu2,
               nfactors=0,distribution="poisson",maxiter=100,maxiter2=0,
               data=borneodata,control=list(trace=1,REPORT=1),method="BFGS",trace=2,estimate=T)


fitinit$log #-8241.699
initmat<-t(matrix(fitinit$coefs[,37],nrow=3)[2:3,])


a<-proc.time()
fitBFGS<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
             convtol=1e-10,maxiter=50,maxiter2=0,init.theta="iterative",init.factor=initmat,
             data=borneodata,control=list(maxit=50,trace=1,REPORT=10),method="BFGS",estimate=T,trace=1)
proc.time()-a

install.packages("maxLik")


a<-proc.time()
fitBFGS<-glmmf2(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
               convtol=1e-10,maxiter=50,maxiter2=0,init.theta="iterative",init.factor=initmat,gradient=F,
               data=borneodata,control=list(maxit=50,trace=1,REPORT=1,fnscale=-1),method="BFGS",estimate=T,trace=1)
proc.time()-a


a<-proc.time()
fitMaxLik<-glmmf2(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
               convtol=1e-10,maxiter=50,maxiter2=0,init.theta="iterative",init.factor=initmat,
               data=borneodata,print.level=3,estimate=T,trace=1,method="NR")
proc.time()-a


#n>30:
#7.77 sec
#1346.198495 


a<-proc.time()
parfit<-parallel_glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
                       convtol=1e-6,maxiter=25,maxiter2=10,init.theta=fitinit$linear.predictor,
                       data=borneodata,control=list(maxit=100,trace=1,REPORT=1),method="BFGS",init.factor=initmat,
                       estimate=T,trace=1,group_size=4,threads=8,samples=8)
proc.time()-a


####



a<-proc.time()
fit<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
            convtol=1e-8,maxiter=50,maxiter2=10,init.theta=fitinit$linear.predictor,
            data=borneodata,control=list(maxit=10000,trace=1),method="Nelder-Mead",init.factor=initmat,estimate=T,trace=1)
proc.time()-a


a<-proc.time()
fit<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
           convtol=1e-8,maxiter=50,maxiter2=10,init.theta=fitinit$linear.predictor,
           data=borneodata,control=list(maxit=1000),method="BFGS",init.factor=initmat,estimate=T,trace=1)
fit$log #-7957.194

microbenchmark(
fit2<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
            convtol=1e-8,maxiter=50,maxiter2=10,init.theta=fitinit$linear.predictor,
            data=borneodata,control=list(maxit=3,trace=1,REPORT=1),method="BFGS",init.factor=initmat,estimate=T,trace=1),times=1)
#58.18007 seconds

split(1:96, cut(1:96, 96%/%12))

a<-proc.time()
fit<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="negative binomial",
           convtol=1e-8,maxiter=50,maxiter2=10,init.theta=fitinit$linear.predictor,
           data=borneodata,control=list(trace=1,REPORT=1,maxit=1000),method="BFGS",init.factor=initmat,estimate=T,trace=1)
fit$log #-7957.194

a<-proc.time()
parfit<-parallel_glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
               convtol=1e-6,maxiter=25,maxiter2=10,init.theta=fitinit$linear.predictor,
               data=borneodata,control=list(maxit=100,trace=1,REPORT=1),method="BFGS",init.factor=initmat,estimate=T,trace=1,group_size=6,threads=8,samples=16,par.maxit=5)
proc.time()-a

parfit10noX<-parfit
save(parfit10noX,file="parfit10noX.rda")

fit<-parfit10noX
fit$log #-1346.198
fit$coefs
factors<-fit$coefs[1:2,]
plot(x=factors[1,],y=factors[2,],pch=NA)
text(x=factors[1,],y=factors[2,], labels=as.character(1:37))

fitinit10noX<-fitinit
fit10noX<-fit
save(fitinit10noX,fit10nox,file="fitBorneo10.rda")

factors<-fit_glmmf2$coefs[1:2,]
dat = read.csv("born.csv")
X = dat[,1:4] #remove id column
plot(x=factors[1,],y=factors[2,],type="n",main="LV + intercept",xlab="LV1",ylab="LV2")
points(x=factors[1,],y=factors[2,],col=as.numeric(X[,2]),pch=as.numeric(X[,3]))
legend("topleft", c("P","L89","L93","L","M","U"),col=c(3,1,2,1,1,1),pch=c(15,15,15,1,2,3))



set.seed(12345)
randu1<-rnorm(37)
randu2<-rnorm(37)
fitinit<-glmmf(group="species",response="abundance",distinct=~Logging+Slope+randu1+randu2,
               nfactors=0,distribution="poisson",maxiter=100,
               data=borneodata,control=list(trace=1,REPORT=1),method="BFGS",trace=2)


glm(abundance~Logging+Slope+randu1+randu2,family=poisson(),data=borneodata[borneodata$species=="Irena_puella",])
fitinit$log #-1433.71
initmat<-t(matrix(fitinit$coefs[,37],nrow=7)[6:7,])



fit<-glmmf(group="species",response="abundance",distinct=~Logging+Slope,nfactors=2,distribution="poisson",
           convtol=1e-8,maxiter=25,maxiter2=10,init.theta=fitinit$linear.predictor,
           data=borneodata,control=list(trace=1,REPORT=1,maxit=10),method="BFGS",init.factor=initmat,estimate=F,trace=2)
fit$log 

(fit$log-fit2$log)/(abs(fit$log)+0.1)
max(abs(fit2$lin-fit$lin))
fit2<-fit
save(fit2,file="fit2.rda")

all.equal(fit$model,fit2$model)
all.equal(fit,fit2)
fit2<-fit
model2<-fit2$model
rnorm(10)
names(fit)
fit$lin
factors<-fit$coefs[1:2,]
plot(x=factors[1,],y=factors[2,],pch=NA)
text(x=factors[1,],y=factors[2,], labels=as.character(1:28))


#####################################################


fitinit<-glmm(response.var="abundance",grouping.var="species",
              fixed=~Logging+Slope+randu1+randu2,
              latent.factors=0,data=borneodata,distribution="poisson")

initmat<-t(matrix(fitinit$fixed$coef,ncol=13)[6:7,])

init.factors<-initmat[lower.tri(initmat,TRUE)]

a<-proc.time()
fit<-glmm(response.var="abundance",grouping.var="species",fixed=~Logging+Slope,
          latent.factors=2,init.factors=init.factors,data=borneodata,distribution="poisson",print_level=1,maxeval=1000000)

proc.time()-a

initmat<-t(tail(matrix(out$coef,ncol=ncol(Y)),2))

init.factors<-initmat[lower.tri(initmat,TRUE)]

Y = dat[,-(1:4)]
X = dat[,2:4] #remove id column
#remove rare species:
Y = as.matrix(Y[,apply(Y>0,2,sum)>10])

Xdat<-do.call(rbind, replicate(ncol(Y), X, simplify=FALSE)) # 5069  3
species<-rep(colnames(Y),each=nrow(Y))
borneodata<-data.frame(c(Y),Xdat,species)
names(borneodata)<-c("abundance",colnames(Xdat),"species")


a<-proc.time()
fit<-glmm(response.var="abundance",grouping.var="species",fixed=~Logging+Slope,
          latent.factors=2,init.factors=init.factors,data=borneodata,distribution="poisson",print_level=1,maxeval=1000000)



fit$fixed$coef
fit$factors$u
save(fit,file="fit30.rda")
plot(unclass(fit$factors$u),pch=NA)
text(unclass(fit$factors$u), labels=as.character(1:37))

plot(unclass(fit$factors$u),type="n",main="LV + main effects",xlab="LV1",ylab="LV2")
points(unclass(fit$factors$u),col=as.numeric(dat[,2]),pch=as.numeric(dat[,3]))
legend("bottomright", c("P","L89","L93","L","M","U"),col=c(3,1,2,1,1,1),pch=c(15,15,15,1,2,3))

####


set.seed(1)
randu1<-rnorm(37)
randu2<-rnorm(37)

fitinit<-glmm(response.var="abundance",grouping.var="species",
              fixed=~randu1+randu2,
              latent.factors=0,data=borneodata,distribution="poisson")

initmat<-t(matrix(fitinit$fixed$coef,ncol=13)[2:3,])

init.factors<-initmat[lower.tri(initmat,TRUE)]

a<-proc.time()
fit<-glmm(response.var="abundance",grouping.var="species",fixed=~1,
          latent.factors=2,init.factors=init.factors,data=borneodata,distribution="poisson",print_level=1,maxeval=1000000)

proc.time()-a


save(fit,file="fit30noX.rda")
plot(unclass(fit$factors$u),type="n",main="LV + intercept",xlab="LV1",ylab="LV2")
points(unclass(fit$factors$u),col=as.numeric(dat[,2]),pch=as.numeric(dat[,3]))
legend("bottomright", c("P","L89","L93","L","M","U"),col=c(3,1,2,1,1,1),pch=c(15,15,15,1,2,3))

####
dat = read.csv("born.csv")
Y = dat[,-(1:4)]
X = dat[,1:4]

XNew <- model.matrix(~Logging+Slope,data=X)
XNew <- XNew[,-1]

#remove rare species:
Y = as.matrix(Y[,apply(Y>0,2,sum)>3])



## Starting values by fitting marginal glm's
start.values <- function(data, X, family, trial.size = NULL, num.lv = 2,seed=1) {
  N <- nrow(data)
  p <- ncol(data)
  data <- as.matrix(data)
  set.seed(seed)
  index <- rmvnorm(N,rep(0,num.lv))
  if(!is.null(X)) fit.mva <- manyglm(data~X+index,family=family, K = trial.size)  
  if(is.null(X)) fit.mva <- manyglm(data~index,family=family, K = trial.size)  
  params <- t(fit.mva$coef)
  if(family == "negative.binomial") { phi <- fit.mva$phi } else { phi <- rep(1,p) }
  
  list(params=params,index=index,phi=phi)
}

library(mvtnorm)
library(mvabund)
sv<-start.values(Y,XNew,"poisson")

plot(ft_main$lvs,type="n",main="LV + main effects",xlab="LV1",ylab="LV2")
points(ft_main$lvs,col=as.numeric(X[,2]),pch=as.numeric(X[,3]))
legend("bottomright", c("P","L89","L93","L","M","U"),col=c(3,1,2,1,1,1),pch=c(15,15,15,1,2,3))
