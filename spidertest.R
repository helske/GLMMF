# Ei toimi, tarvitsee tasoituksen faktoreita varten, nehän riippuu ajasta!
# Voiko tasoitusta pelkistää niin että tehdään se vain faktoritiloille?
# jossain vaiheessa betalle priori, oletus eksakti diffuusi, mutta voi laittaa myös var<inf ja E(beta) jotain
# katso rosenin artikkelia vielä tarkemmin, pitäisikö rakentaa suoraan rinnakkaislaskennallinen versio?
# erikseen funktiot mallille jossa ei latentteja faktoreita ja erikseen niiden kanssa.

load("../gitrepos/GLMMF/kmodel.rda")

load("../StateSpaceGLMM/spiderdata.rda")
library(GLMMF)

set.seed(1)
randu1<-rnorm(max(spiderdata[,"site"]))
randu2<-rnorm(max(spiderdata[,"site"]))

fitinit<-glmmf(group="species",response="abundance",distinct=~randu1+randu2,nfactors=0,distribution="poisson",
           data=spiderdata,control=list(trace=1,REPORT=1),method="BFGS")

fitinit$log #-3307.467

initmat<-t(matrix(fitinit$coefs[,28],nrow=3)[2:3,])



fit<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
           data=spiderdata,control=list(trace=1,REPORT=1,maxit=1000),method="BFGS",init.factor=initmat,estimate=TRUE,trace=1)

factors<-fit$coefs[1:2,]
plot(x=factors[1,],y=factors[2,],pch=NA)
text(x=factors[1,],y=factors[2,], labels=as.character(1:28))

########

approxSSM(kmodel)$iter

init.factors<-initmat[lower.tri(initmat,TRUE)]
fit2<-glmm(response.var="abundance",grouping.var="species",fixed=~1,
          latent.factors=2,init.factors=init.factors,data=spiderdata,distribution="poisson",print_level=2,maxeval=100000)



fit<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",maxiter=0,
           data=spiderdata,control=list(trace=1,REPORT=1),method="BFGS",init.factor=initmat,estimate=FALSE)

debug(glmmf)
x<-matrix(0,12,2)
x[,1]<-pars[1:12]
x[2:12,2]<-pars[13:23]
kmodel$Z[,13:14,]<-x
logLik(kmodel)
fit<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
           data=spiderdata,control=list(trace=1,REPORT=1),method="Nelder-Mead",init.factor=x,estimate=TRUE)


suppressWarnings(deviance(KFS(approxSSM(kmodel,maxiter=2))))

dev<-rep(0,15)
for(i in 0:14)
dev[i+1]<-suppressWarnings(deviance(KFS(kmodel,maxiter=i)))

lla<-rep(0,15)
for(i in 0:14)
  lla[i+1]<-suppressWarnings(logLik(approxSSM(kmodel,maxiter=i)))

ll<-rep(0,15)
for(i in 0:14)
  ll[i+1]<-suppressWarnings(logLik(kmodel,maxiter=i))

ts.plot(dev)
ts.plot(ll)
ts.plot(dev-ll)
diff(dev)
diff(ll)
suppressWarnings(deviance(KFS(kmodel,maxiter=2)))
load("kmodel.rda")
load("gmodel.rda")
model<-gmodel

load("pars.rda")
initmat<-matrix(0,12,2)
initmat[-13]<-pars

init.factors2<-initmat[lower.tri(initmat,TRUE)]

kmodel$Z[,13:14,]<-initmat
model$Z[1:2,,]<-t(initmat)



theta <- KFAS:::init_theta(model$y, model$u, rep(model$distribution,12)) 
maxiter<-15
#maxiter<-maxiter+1
res<-GLMMF:::expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                     rep(2,12), model$tol, maxiter, 1e-10, theta, model$Zind, model$nfactors)


dev[maxiter]
lla[maxiter]
res$log+KFAS:::scaling(model$y,model$u,rep(model$distribution,12),res$theta)
ll[maxiter]
res$conv
app<-approxSSM(kmodel,maxiter=0)
all.equal(res$theta,unclass(app$theta),check.attributes=FALSE)
o<-KFS(approxSSM(kmodel,maxiter=0),simp=F)
o$r[,28]
app$theta[1:5,1:5]
res$theta[1:5,1:5]

init.factors2<-initmat[lower.tri(initmat,TRUE)]

fit<-glmm(response.var="abundance",grouping.var="species",fixed=~1,
          latent.factors=2,init.factors=init.factors2,data=spiderdata,distribution="poisson",print_level=3,maxeval=0)


all.equal(fitinit$coefs[,28],out$alpha[28,],check.attributes=FALSE)
library(KFAS)
kmodel<-SSModel(model$y~randu1+randu2,dist="poisson")
out<-KFS(kmodel,smoothing=c("state","mean","signal"))

out<-approxSSM(kmodel,maxiter=100)
out$theta-fitinit$fitted
kmodel$y-model$y
logLik(kmodel)
fitinit$log
init.factors<-initmat[lower.tri(initmat,TRUE)]

fit<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
           data=spiderdata,control=list(trace=1,REPORT=1),method="BFGS")


data(butterfly)
fit.SSglmm<-glmm(response.var="Colias", grouping.var="site", fixed=~habitat + building + urbanveg,random=~1,
                 data=butterfly,distribution="poisson")

fit.SSglmm2<-glmm(response.var="Colias", grouping.var="site", fixed=~habitat + building + urbanveg,random=~1,
                  latent.factors=2,data=butterfly,distribution="poisson",print_level=1)

data(spider)


Y<-as.matrix(spider$abund)
X<-spider$x
X<-do.call(rbind, replicate(ncol(Y), X, simplify=FALSE))
species<-rep(colnames(spider$abund),each=nrow(Y))
spiderdata<-data.frame(c(Y),X,species,site=rep(1:nrow(Y),ncol(Y)))
names(spiderdata)<-c("abundance",colnames(X),"species","site")
save(spiderdata,file="spiderdata.rda")
debug(glmm)
load("spiderdata.rda")


set.seed(2)
randu1<-rnorm(max(spiderdata[,"site"]))
randu2<-rnorm(max(spiderdata[,"site"]))

fitinit<-glmm(response.var="abundance",grouping.var="species",
              fixed=~1+randu1+randu2,
              latent.factors=0,data=spiderdata,distribution="poisson")
initmat<-t(matrix(fitinit$fixed$coef,nrow=3)[2:3,])

init.factors<-initmat[lower.tri(initmat,TRUE)]

fit<-glmm(response.var="abundance",grouping.var="species",fixed=~1,
          latent.factors=2,init.factors=init.factors,data=spiderdata,distribution="poisson",print_level=1,maxeval=1000000)

fit$fixed$coef
fit$factors$u

plot(x=as.numeric(fit$factors$u[,1]),y=as.numeric(fit$factors$u[,2]),pch=NA)
text(x=as.numeric(fit$factors$u[,1]),y=as.numeric(fit$factors$u[,2]), labels=as.character(1:28))

fit1<-fit
fitinit1<-fitinit
####

datainit<-cbind(spiderdata,randu1,randu2)
i=5
summary(glm(abundance~soil.dry+bare.sand+fallen.leaves+randu1,
            data=datainit[datainit[,"species"]==colnames(spider$abund)[i],],family=poisson()))

fit<-manyglm(Y~as.matrix(datainit[1:28,c(2:4,10)]),family="poisson",show.warning=TRUE,
             tol=1e-8)
fit$iter
fit$coef
fit$fitted ## T??ll? on jotain rajoitteita ja siksi ei herjaa mit??n toisin kuin glm!!!

latent.factors<-2
set.seed(123)
randu1<-rnorm(max(spiderdata[,"site"]))
randu2<-rnorm(max(spiderdata[,"site"]))

datainit<-cbind(spiderdata,randu1,randu2)
i=5
summary(glm(abundance~soil.dry+bare.sand+fallen.leaves+moss+herb.layer+reflection + randu1+randu2,
            data=datainit[datainit[,"species"]==colnames(spider$abund)[i],],family=poisson()))

fit<-manyglm(Y~as.matrix(datainit[1:28,c(2:7,10,11)]),family="poisson")
fit$coef

i<-0
i<-i+1
summary(glm(abundance~soil.dry+bare.sand+fallen.leaves+moss+herb.layer+reflection + randu1+randu2,
    data=spiderdata[spiderdata[,"species"]==colnames(spider$abund)[i],],family=poisson()))


#i=5

datainit<-cbind(spiderdata,randu1,randu2)
i=5
colnames(datainit)
model<-SSModel(abundance~soil.dry+bare.sand+fallen.leaves+moss+herb.layer+reflection + randu1+randu2,
               data=spiderdata[spiderdata[,"species"]==colnames(spider$abund)[i],],distribution="poisson")
app<-approxSSM(model,maxiter=1000)
out<-KFS(app)
out$d
out$V[,,28]-out$V[,,27]


i<-i+1

model<-SSModel(abundance~soil.dry+moss+herb.layer+reflection + randu1+randu2,
               data=spiderdata[spiderdata[,"species"]==colnames(spider$abund)[i],],distribution="poisson")

app<-approxSSM(model,maxiter=100)
out<-KFS(app)
out$d
sqrt(diag(out$V[,,28]))
KFS(SSModel(abundance~soil.dry+moss+herb.layer+reflection + randu1+randu2,
    data=spiderdata[spiderdata[,"species"]==colnames(spider$abund)[i],],distribution="poisson",tol=1e-8),maxiter=1000)




set.seed(12345)
randu1<-rnorm(max(spiderdata[,"site"]))
randu2<-rnorm(max(spiderdata[,"site"]))

fitinit<-glmm(response.var="abundance",grouping.var="species",
              fixed=~soil.dry+moss+herb.layer+reflection + randu1+randu2,
          latent.factors=0,data=spiderdata,distribution="poisson",print_level=3,maxeval=1000000)

fit<-glmm(response.var="abundance",grouping.var="species",fixed=~soil.dry+moss+herb.layer+reflection,
          latent.factors=latent.factors,data=spiderdata,distribution="poisson",print_level=3,maxeval=1000000)


names(fit)
fit$fixed$coef
fit$fixed
###
fit2<-fit

fit<-glmm(response.var="abundance",grouping.var="species",fixed=~soil.dry,
          latent.factors=0,data=spiderdata,distribution="poisson",print_level=3)
out<-KFS(fit$model,smoothing=c("state","mean"))

app<-approxSSM(fit$model)
outapp<-KFS(app)
outapp$Finf
dim(fit$model$y)
