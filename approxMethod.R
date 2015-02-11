library(GLMMF)
library(StateSpaceGLMM)
library(nloptr)
library(KFAS)
library(numDeriv)
dat = read.csv("../GLMMF/born.csv")
Y = dat[,-(1:4)]
X = dat[,2:4] #remove id column
#remove rare species:
Y = as.matrix(Y[,apply(Y>0,2,sum)>30])

Xdat<-do.call(rbind, replicate(ncol(Y), X, simplify=FALSE)) # 5069  3
species<-rep(colnames(Y),each=nrow(Y))
borneodata<-data.frame(c(Y),Xdat,species)
names(borneodata)<-c("abundance",colnames(Xdat),"species")

set.seed(12345)
randu1<-rnorm(37)
randu2<-rnorm(37)

fitinit_glmmf<-glmmf(group="species",response="abundance",distinct=~randu1+randu2,
                     nfactors=0,distribution="poisson",maxiter=100,
                     data=borneodata,control=list(trace=1,REPORT=1),method="BFGS",trace=2,estimate=F)

fitinit_glmmf$log #-1561.845
initmat_glmmf<-t(matrix(fitinit_glmmf$coefs[,37],nrow=3)[2:3,])

a<-proc.time()
fit_glmmf2<-glmmf2(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
                 convtol=1e-8,maxiter=100,maxiter2=5,init.factor=initmat_glmmf,
                 data=borneodata, estimate=T,trace=1,epsilon=1e-8,maxapp=100,gradient=TRUE,
                 check.analyticals=FALSE,print.level=0,typsize=rep(0,sum(lower.tri(initmat_glmmf,TRUE))),
                 fscale=7000,iterlim=50)
proc.time()-a
# -7127.43818564259
# 1301.23 sec
# 7 maxapp-iteraatiota
#

a<-proc.time()
fit_glmmf3<-glmmf2(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
                   convtol=1e-8,maxiter=100,maxiter2=5,init.factor=initmat_glmmf,
                   data=borneodata, estimate=T,trace=1,epsilon=1e-8,maxapp=100,gradient=TRUE,
                   check.analyticals=FALSE,print.level=0,typsize=rep(0,sum(lower.tri(initmat_glmmf,TRUE))),
                   fscale=7000,iterlim=500)
proc.time()-a
# -7127.43818564259
# 1791.86 
# "ApproxMaximization iter 6 log-likelihood -7127.43818564259"

load("fit10approx.rda")

factors<-fit10a$coefs[1:2,]
dat = read.csv("born.csv")
X = dat[,1:4] #remove id column
plot(x=factors[1,],y=factors[2,],type="n",main="LV + intercept",xlab="LV1",ylab="LV2")
points(x=factors[1,],y=factors[2,],col=as.numeric(X[,2]),pch=as.numeric(X[,3]))
legend("bottomright", c("P","L89","L93","L","M","U"),col=c(3,1,2,1,1,1),pch=c(15,15,15,1,2,3))

t(fit10a$model$Z[1:2,,1])

#fit_glmmf2$logLik
# >10: -7957.194

# >20:
# "Iter 4 log-likelihood 4628.61211358786"
# > proc.time()-a
# user  system elapsed 
# 128.13    0.32  128.81 

library(microbenchmark)

microbenchmark(
fit_glmmf2<-glmmf2(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
                   convtol=1e-8,maxiter=100,maxiter2=5,init.factor=initmat_glmmf,
                   data=borneodata, estimate=F,trace=0,epsilon=1e-5,print.level=0,gradient=TRUE,check.analyticals=FALSE),

fit_glmmf2<-glmmf2(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
                   convtol=1e-8,maxiter=100,maxiter2=5,init.factor=initmat_glmmf,
                   data=borneodata, estimate=T,trace=0,epsilon=1e-5,print.level=0,gradient=TRUE,check.analyticals=TRUE),

fit_glmmf2<-glmmf2(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
                   convtol=1e-8,maxiter=100,maxiter2=5,init.factor=initmat_glmmf,
                   data=borneodata, estimate=T,trace=0,epsilon=1e-5,print.level=0,gradient=FALSE),times=5)
# min         lq       mean     median         uq        max neval cld
# 686.1772   691.2514   692.6971   692.0845   692.8338   701.1388     5 a  
# 1119.0128  1121.8169  1125.9374  1123.4147  1129.3156  1136.1267     5  b 
# 11880.2862 11916.1517 11939.8279 11947.8354 11974.2710 11980.5952     5   c
fit_glmmf2$logLik #-1346.198


GLMMF:::expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                2, model$tol, 100, 0, 1e-5, theta, model$Zind, model$nfactors,trace,1)$logLik
