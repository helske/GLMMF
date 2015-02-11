library(KFAS)
load("spider.rda")
set.seed(1)
randu1<-rnorm(28)
randu2<-c(0,rnorm(27))

initmodel<-SSModel(ts(spider$abund) ~ randu1+randu2
                   ,data=data.frame(spider$x,randu1=randu1,randu2=randu2),distribution="poisson")

initout<-KFS(initmodel)
Z<-matrix(coef(initout,start=1,end=1),ncol=3,byrow=TRUE)[,2:3]
Z[upper.tri(Z)]<-0

model<-SSModel(ts(spider$abund) ~
                 SSMcustom(Z=Z,T=diag(0,2),
                           R=diag(2),Q=diag(2),P1=diag(2),
                           P1inf=diag(0,2),a1=matrix(0,2)),
               ,data=data.frame(spider$x),distribution="poisson")


trace(KFAS:::logLik.SSModel, print=FALSE,exit=
        quote({
          if(all(object$distribution == "gaussian")){
            .countG <<- .countG + 1
          } else {
            if(out$info == 1) .countNaN <<- .countNaN+1
            if(out$info == 2) .countIter <<- .countIter+1
            .countNG <<- .countNG + out$maxiter
          }         
        }))
trace(KFAS:::KFS, print=FALSE,exit=
        quote({
            .countG <<- .countG + 1
          } ))

likfn<-function(pars,model,estimate=TRUE){
  Z<-matrix(0,12,2)
  Z[lower.tri(Z,TRUE)]<-pars
  model$Z[,13:14,]<-Z   
  -KFAS:::logLik.SSModel(model,check.model=FALSE)    
}

.countG<-.countNG<-.countNaN<-.countIter<-0
fit1<-optim(par=model$Z[,13:14,1][lower.tri(matrix(0,12,2),TRUE)],
            fn=likfn, model=model,method='BFGS',control=list(REPORT=1,trace=1,maxit=10000))
standardmethod<-c(.countG,.countNG,.countNaN,.countIter) # 0 218990      0      2
model1<-model
Z1<-matrix(0,12,2)
Z1[lower.tri(Z1,TRUE)]<-fit1$par
model1$Z[,13:14,]<-Z1
Z1

out1<-KFS(model1)
lv1<-unclass(coef(out1)[,13:14])
plot(x=lv1[,1],y=lv1[,2],pch=NA)
text(x=lv1[,1],y=lv1[,2],labels=as.character(1:28))

grad<-function(pars,model){
  Z<-matrix(0,12,2)
  Z[lower.tri(Z,TRUE)]<-pars
  model$Z[,13:14,]<-Z   
  p<-attr(model,"p")
  n<-attr(model,"n")
  out<-KFS(model,filtering="none",smoothing="state")
  gradient<-matrix(0,p,2)
  
  for(j in 1:p){
    for(k in 1:2){
      if(j>=k){
        for(t in 1:n){
          gradient[j,k] <- gradient[j,k] + 
            (model$y[t,j]*out$alpha[t,12+k]-(out$V[12+k,,t]+out$alpha[t,12+k]%*%t(out$alpha[t,]))%*%model$Z[j,,t])/model$H[j,j,t]
        }
      }
    }
  }
   -gradient[lower.tri(gradient,TRUE)]
}



.countG<-.countNG<-.countNaN<-.countIter<-0
model2<-model
for(i in 1:10){
  amodel<-approxSSM(model2)
  .countNG<-.countNG+amodel$iterations
  fita<-optim(par=model2$Z[,13:14,1][lower.tri(matrix(0,12,2),TRUE)],
              fn=likfn, gr=grad,model=amodel,method='BFGS',control=list(trace=1,maxit=10))
  Z2a<-matrix(0,12,2)
  Z2a[lower.tri(Z2a,TRUE)]<-fita$par
  model2$Z[,13:14,]<-Z2a
}
fita<-optim(par=model2$Z[,13:14,1][lower.tri(matrix(0,12,2),TRUE)],
            fn=likfn, gr=grad, model=amodel,method='BFGS',control=list(REPORT=1,trace=1,maxit=10000))
fit2<-optim(par=fita$p,
            fn=likfn, model=model2,method='BFGS',control=list(REPORT=1,trace=1,maxit=10000))
newmethod<-c(.countG,.countNG,.countNaN,.countIter) #312 21669     0     0
Z2<-matrix(0,12,2)
Z2[lower.tri(Z2,TRUE)]<-fit2$par
model2$Z[,13:14,]<-Z2
Z2


gradloglik<-function(pars,model){
  p<-attr(model,"p")
  n<-attr(model,"n")
  out<-KFS(model,filtering="state",smoothing="state")
  gradient<-matrix(0,p,2)
  
  for(j in 1:p){
    for(k in 1:2){
      if(j>=k){
        for(t in 1:n){
          gradient[j,k] <- gradient[j,k] + 
            (model$y[t,j]*out$alpha[t,12+k]-(out$V[12+k,,t]+out$alpha[t,12+k]%*%t(out$alpha[t,]))%*%model$Z[j,,t])/model$H[j,j,t]
        }
      }
    }
  }
  ll<- -out$logLik
  attr(ll,"gradient")<- -gradient[lower.tri(gradient,TRUE)]
  ll
}


