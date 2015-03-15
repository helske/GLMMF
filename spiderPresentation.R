library(KFAS)
load("spider.rda")
set.seed(1)
randu1<-rnorm(28)
randu2<-c(0,rnorm(27))

initmodel<-SSModel(ts(spider$abund) ~ randu1+randu2
                   ,data=data.frame(spider$x,randu1=randu1,randu2=randu2),distribution="poisson")
logLik(initmodel) # -3436.87
initout<-KFS(initmodel)
Z<-matrix(coef(initout,start=1,end=1),ncol=3,byrow=TRUE)[,2:3]
Z[upper.tri(Z)]<-0

model<-SSModel(ts(spider$abund) ~
                 SSMcustom(Z=Z,T=diag(0,2),
                           R=diag(2),Q=diag(2),P1=diag(2),
                           P1inf=diag(0,2),a1=matrix(0,2)),
               ,data=data.frame(spider$x),distribution="poisson")
logLik(model) #-2256.94
# Count the number of gaussian and non-gaussian filtering (and smoothing) were performed
trace(KFAS:::logLik.SSModel, print=FALSE,exit=
        quote({
          if(all(object$distribution == "gaussian")){
            .countG <<- .countG + 1
          } else {
            if(out$info %in% (1:2) ) .countNaN <<- .countNaN+1
            if(out$info == 3) .countIter <<- .countIter+1
            .countNG <<- .countNG + out$maxiter
          }
        }))


likfn<-function(pars,model,estimate=TRUE){
  Z<-matrix(0,12,2)
  Z[lower.tri(Z,TRUE)]<-pars
  model$Z[,13:14,]<-Z
  -KFAS:::logLik.SSModel(model,check.model=FALSE)
}

.countG<-.countNG<-.countNaN<-.countIter<-0
fit1<-optim(par=model$Z[,13:14,1][lower.tri(matrix(0,12,2),TRUE)],
            fn=likfn, model=model,method='BFGS',control=list(REPORT=1,trace=1,maxit=1000))
fit1$val#857.2482
standardway<-c(.countG,.countNG,.countNaN,.countIter) #      0 218990  0      2
model1<-model
Z1<-matrix(0,12,2)
Z1[lower.tri(Z1,TRUE)]<-fit1$par
model1$Z[,13:14,]<-Z1
Z1

out1<-KFS(model1)
lv1<-unclass(coef(out1)[,13:14])
plot(x=lv1[,1],y=lv1[,2],pch=NA)
text(x=lv1[,1],y=lv1[,2],labels=as.character(1:28),col=2)


.countG<-.countNG<-.countNaN<-.countIter<-0
model2<-model
for(i in 1:10){
  amodel<-approxSSM(model2)
  .countNG<-.countNG+amodel$iterations
  fita<-optim(par=model2$Z[,13:14,1][lower.tri(matrix(0,12,2),TRUE)],
              fn=likfn, model=amodel,method='BFGS',control=list(trace=1,maxit=10))
  Z2a<-matrix(0,12,2)
  Z2a[lower.tri(Z2a,TRUE)]<-fita$par
  model2$Z[,13:14,]<-Z2a
}
fita<-optim(par=model2$Z[,13:14,1][lower.tri(matrix(0,12,2),TRUE)],
            fn=likfn, model=amodel,method='BFGS',control=list(REPORT=1,trace=1,maxit=10000))
fit2<-optim(par=fita$p,
            fn=likfn, model=model2,method='BFGS',control=list(REPORT=1,trace=1,maxit=10000))
newmethod<-c(.countG,.countNG,.countNaN,.countIter) #9199 21844     0     0
Z2<-matrix(0,12,2)
Z2[lower.tri(Z2,TRUE)]<-fit2$par
model2$Z[,13:14,]<-Z2
Z2

out2<-KFS(model2)
lv2<-unclass(coef(out2)[,13:14])
plot(x=lv2[,1],y=lv2[,2],pch=NA)
text(x=lv2[,1],y=lv2[,2],labels=as.character(1:28))
text(x=lv1[,1],y=lv1[,2],labels=as.character(1:28),col=2)
lv1-lv2
par(mfrow=c(2,1))
logLik(model1)
logLik(model2)

######## Negative binomial


initmodelNB<-SSModel(ts(spider$abund) ~ randu1+randu2
                     ,data=data.frame(cbind(spider$x,randu1=randu1,randu2=randu2)),distribution="negative binomial")
initoutNB<-KFS(initmodelNB)

ZNB<-matrix(coef(initoutNB,start=1,end=1),ncol=3,byrow=TRUE)[,2:3]
ZNB[upper.tri(ZNB)]<-0


modelNB<-SSModel(ts(spider$abund) ~
                   SSMcustom(Z=Z,T=diag(0,2),
                             R=diag(2),Q=diag(2),P1=diag(2),
                             P1inf=diag(0,2),a1=matrix(0,2)),
                 ,data=data.frame(spider$x),distribution="negative binomial")

outNB<-KFS(modelNB)

likfnNB<-function(pars,model,estimate=TRUE){
  Z<-matrix(0,12,2)
  Z[lower.tri(Z,TRUE)]<-pars[1:23]
  model$Z[,13:14,]<-Z
  model$u[]<-rep(exp(pars[24:35]),each=28)
  -KFAS:::logLik.SSModel(model,check.model=FALSE)
}

a<-proc.time()
.countG<-.countNG<-.countNaN<-.countIter<-0
fit1NB<-optim(par=c(model$Z[,13:14,1][lower.tri(matrix(0,12,2),TRUE)],rep(0,12)),
              fn=likfnNB, model=modelNB,method='BFGS',control=list(maxit=10000))
standardwayNB<-c(.countG,.countNG,.countNaN,.countIter) # 0 64015     0   145
proc.time()-a

model1NB<-modelNB
Z1NB<-matrix(0,12,2)
Z1NB[lower.tri(Z1,TRUE)]<-fit1NB$par[1:23]
model1NB$Z[,13:14,]<-Z1NB
Z1NB
model1NB$u[]<-rep(exp(fit1NB$par[24:35]),each=28)
out1NB<-KFS(model1NB)

lv1NB<-unclass(coef(out1NB)[,13:14])
plot(x=lv1NB[,1],y=lv1NB[,2],pch=NA,ylim=range(c(lv1NB[,2],lv1[,2])),xlim=range(c(lv1NB[,1],-lv1[,1])))
text(x=lv1NB[,1],y=lv1NB[,2],labels=as.character(1:28))
text(x=-lv1[,1],y=lv1[,2],labels=as.character(1:28),col=3) #minus sign for rotation

logLik(model1NB)
logLik(model1)

.countG<-.countNG<-.countNaN<-.countIter<-0
model2NB<-modelNB
for(i in 1:2){
  amodelNB<-approxSSM(model2NB)
  .countNG<-.countNG+amodelNB$iterations
  fitaNB<-optim(par=c(model2NB$Z[,13:14,1][lower.tri(matrix(0,12,2),TRUE)]),
                fn=likfn, model=amodelNB,method='BFGS',control=list(maxit=10))
  Z2aNB<-matrix(0,12,2)
  Z2aNB[lower.tri(Z2a,TRUE)]<-fitaNB$par
  model2NB$Z[,13:14,]<-Z2aNB
}

fitaNB<-optim(par=c(model2NB$Z[,13:14,1][lower.tri(matrix(0,12,2),TRUE)]),
              fn=likfn, model=amodelNB,method='BFGS',control=list(maxit=10000))
fit2NB<-optim(par=c(fitaNB$p,rep(0,12)),
              fn=likfnNB, model=model2NB,method='BFGS',control=list(maxit=10000))
newmethodNB<-c(.countG,.countNG,.countNaN,.countIter) # 2125 39546     0     4
model2NB<-modelNB
Z2NB<-matrix(0,12,2)
Z2NB[lower.tri(Z2NB,TRUE)]<-fit2NB$par[1:23]
model2NB$Z[,13:14,]<-Z2NB
Z2NB
model2NB$u[]<-rep(exp(fit2NB$par[24:35]),each=28)
out2NB<-KFS(model2NB)
logLik(model1NB)
logLik(model2NB)

lv2NB<-unclass(coef(out2NB)[,13:14])
lv1NB<-unclass(coef(out1NB)[,13:14])
plot(x=lv1NB[,1],y=lv1NB[,2],pch=NA,ylim=range(c(lv1NB[,2],lv2NB[,2])),xlim=range(c(lv1NB[,1],lv2NB[,1])))
text(x=lv1NB[,1],y=lv1NB[,2],labels=as.character(1:28))
text(x=lv2NB[,1],y=lv2NB[,2],labels=as.character(1:28),col=2)
out2NBs<-KFS(model2NB,nsim=1000)
lv2NBs<-unclass(coef(out2NBs)[,13:14])
text(x=lv2NBs[,1],y=lv2NBs[,2],labels=as.character(1:28),col=3)


###############
