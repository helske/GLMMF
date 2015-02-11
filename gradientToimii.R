# Poisson model
model<-SSModel(VanKilled~-1+SSMcustom(Z=1,T=1,a1=1,P1=0,P1inf=0,Q=0,R=0)+law+
                 SSMseasonal(period=12,sea.type='dummy',Q=0),
               data=Seatbelts, distribution='poisson')

likfn<-function(psi){
  model$Z[,13,]<-psi
  set.seed(1)
  -logLik(model,tol=1e-15,maxiter=1000)
}

gradf<-function(psi){
  model$Z[,13,]<-psi
  set.seed(1)
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  -sum(out$alpha[,13]*(model$y-exp(out$theta)))
}

#model$P1inf[]<-0
#model$P1inf[1]<-1
grad(likfn,2)-gradf(2)


f<-function(psi){
  model$Z[,13,]<-psi
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  app<-approxSSM(model)
  -sum(dpois(model$y,exp(out$theta),log=TRUE))+sum(dnorm(app$y,out$theta,sqrt(c(apply(app$H,3,diag))),log=TRUE))-logLik(app)
}

grad(likfn,2)
grad(f,2)


fpois<-function(psi){
  model$Z[,13,]<-psi
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  app<-approxSSM(model)
  -sum(dpois(model$y,exp(out$theta),log=TRUE))
}

fnorm<-function(psi){
  model$Z[,13,]<-psi
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  app<-approxSSM(model)
  sum(dnorm(app$y,out$theta,sqrt(c(apply(app$H,3,diag))),log=TRUE))
}

fapp<-function(psi){
  model$Z[,13,]<-psi
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  app<-approxSSM(model)
  -logLik(app)
}


gpois<-function(psi){
  model$Z[,13,]<-psi
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  app<-approxSSM(model)
  -sum(out$alpha[,13]*(model$y-exp(out$theta)))
}


gnorm<-function(psi){
  model$Z[,13,]<-psi
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  app<-approxSSM(model)
  sum(-0.5*(exp(app$theta)-model$y^2*exp(-app$theta)-1)*out$alpha[,13])
  
}


gapp<-function(psi){
  model$Z[,13,]<-psi
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  app<-approxSSM(model)
  g<-0
  for(i in 1:nrow(model$y))
    g <- g + (app$y[i]*out$alpha[i,13]-out$V[13,,i]%*%model$Z[1,,i]-out$alpha[i,13]%*%t(out$alpha[i,])%*%model$Z[1,,i])/app$H[1,1,i]
  g
}

grad(fpois,2)+grad(fnorm,2)+grad(fapp,2)


gnorm(2)
gpois(2)
gapp(2)

grad(likfn,2)
grad(f,2)



fnorm2<-function(psi){
  model$Z[,13,]<-psi
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  app<-approxSSM(model)
  sum(-0.5*(-app$theta+model$y^2*exp(-app$theta)+exp(app$theta)-2*model$y))
}

gnorm<-function(psi){
  model$Z[,13,]<-psi
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  app<-approxSSM(model,tol=1e-15,maxiter=1000)
  sum(-0.5*(exp(app$theta)-model$y^2*exp(-app$theta)-1)*out$alpha[,13])
  
}


gnorm<-function(psi){
  model$Z[,13,]<-psi
  app<-approxSSM(model,tol=1e-15,maxiter=1000)
  out<-KFS(app,smoothing=c("state","signal","disturbance"))  
   - sum(0.5*(1+2*app$theta*exp(-app$theta)-app$theta^2*exp(-app$theta))*out$alpha[,13])  
}
grad(fnorm2,2)
gnorm(2)
grad(fapp,2)
grad(fnorm,2)
all.equal(grad(fnorm2,2),grad(fnorm,2))


gradf<-function(psi){
  model$Z[,13,]<-psi
  app<-approxSSM(model,tol=1e-15,maxiter=1000)
  out<-KFS(app)
  #out<-KFS(model,smoothing=c("state","signal"),maxiter=0) #,convtol=1e-15,maxiter=1000)
  g<-0
  for(i in 1:nrow(model$y))
    g <- g - (out$alpha[i,]*model$y[i]-(out$alpha[i,]*exp(app$theta[i])+0.5*exp(app$theta[i])*(app$theta[i]+2)*model$Z[1,,i]%*%out$V[,,i]))
  g
}
gradf(2)
fit<-optim(f=likfn,gr=gradf,p=2,method="BFGS")
model$Z[,13,]<-fit$par

out<-KFS(model)

ts.plot(cbind(model$y,out$mu),col=1:2)


diag(model$P1)[1:12]<-100000




likfn<-function(psi){
  model$Z[,13,]<-psi
  -logLik(model,maxiter=0)
}






likfn2<-function(psi){
  model$Z[,13,]<-psi
  lik<- -logLik(model,tol=1e-15,maxiter=1000)
  attr(lik,"gradient")<-unname(gradf(psi))
  lik
}




gradf(2)
grad(likfn,2,method.args=list(r=6))

fit2<-optim(f=likfn,gr=gradf,p=2,method="BFGS")

f<-function(x) Vectorize(gradf)(x)-grad(Vectorize(likfn),x)
curve(f(x),from=0,to=3)
curve(grad(Vectorize(likfn),x),from=0,to=3)

root<-uniroot(gradf,interval=c(0,3),tol=1e-15)
gradf(root$root)

f<-function(x){ grad(likfn,x)}
uniroot(f,interval=c(0,3))

likfnP<-function(psi){
  model$Z[1,13,]<-psi
  out<-KFS(model,smoothing=c("signal"),convtol=1e-15,maxiter=1000)
  -sum(dpois(model$y,exp(out$theta),log=TRUE))
}

grad(likfnP,fit2$par)-grad(likfn,fit2$par)
gradf(10)
