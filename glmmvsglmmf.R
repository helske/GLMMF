library(GLMMF)
library(StateSpaceGLMM)
library(nloptr)
library(KFAS)
library(numDeriv)
dat = read.csv("../GLMMF/born.csv")
Y = dat[,-(1:4)]
X = dat[,2:4] #remove id column
#remove rare species:
Y = as.matrix(Y[,apply(Y>0,2,sum)>10])

Xdat<-do.call(rbind, replicate(ncol(Y), X, simplify=FALSE)) # 5069  3
species<-rep(colnames(Y),each=nrow(Y))
borneodata<-data.frame(c(Y),Xdat,species)
names(borneodata)<-c("abundance",colnames(Xdat),"species")
save(borneodata,file="borneodata.rda")

load("borneodata.rda")

set.seed(12345)
randu1<-rnorm(37)
randu2<-rnorm(37)

fitinit_glmm<-glmm(response.var="abundance",grouping.var="species",
              fixed=~randu1+randu2,
              latent.factors=0,data=borneodata,distribution="poisson",maxiter=100)
fitinit_glmm$log #-1561.845
initmat_glmm<-t(matrix(fitinit_glmm$fixed$coef,ncol=13)[2:3,])

fitinit_glmmf<-glmmf(group="species",response="abundance",distinct=~randu1+randu2,
               nfactors=0,distribution="poisson",maxiter=100,
               data=borneodata,control=list(trace=1,REPORT=1),method="BFGS",trace=2,estimate=F)


fitinit_glmmf$log #-1561.845
initmat_glmmf<-t(matrix(fitinit_glmmf$coefs[,37],nrow=3)[2:3,])

all.equal(unname(fitinit_glmm$fixed$coef),fitinit_glmmf$coefs[,37])



fit_glmm<-glmm(response.var="abundance",grouping.var="species",fixed=~1,maxiter=100,convtol=1e-10,
                                latent.factors=2,init.factors=initmat_glmm[lower.tri(initmat_glmm,TRUE)],
               data=borneodata,distribution="poisson",print_level=1,estimate=T,maxeval=100000)
fit_glmm$logLik #-1346.198489

fit_glmm_true<-fit_glmm
fit_glmm_true$logLik

amodel<-approxSSM(fit_glmm$model)

likfn<-function(pars,model,estimate=TRUE){
  Z<-matrix(0,13,2)
  Z[lower.tri(Z,TRUE)]<-pars
  model$Z[,14:15,]<-Z
  if(estimate)
    -logLik(model)
  else model
}

fit<-optim(f=likfn,p=initmat_glmm[lower.tri(initmat_glmm,TRUE)],model=amodel,method="BFGS",control=list(trace=1,REPORT=1))

model<-fit_glmm$model
Z<-matrix(0,13,2)
Z[lower.tri(Z,TRUE)]<-fit$p
model$Z[,14:15,]<-Z
logLik(model)
fit<-optim(f=likfn,p=Z[lower.tri(Z,TRUE)],model=approxSSM(model),method="BFGS",control=list(trace=1,REPORT=1))

Z<-matrix(0,13,2)
Z[lower.tri(Z,TRUE)]<-fit$p
model$Z[,14:15,]<-Z
logLik(model)

fit_glmmf<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
                 convtol=1e-8,maxiter=100,maxiter2=5,init.theta="iterative",init.factor=initmat_glmmf,
                 data=borneodata,control=list(trace=1,REPORT=1,maxit=100),method="BFGS",
                 estimate=T,trace=1)
fit_glmmf$logLik #-1346.198

model<-fit_glmm$model

likfn1<-function(pars){
model$Z[1,14,]<-pars
logLik(model,convtol=1e-15,maxiter=500)
}

grad(likfn1,1)



gpois<-function(psi){
  model$Z[1,14,]<-psi
  out<-KFS(model,smoothing=c("state","signal"),convtol=1e-15,maxiter=1000)
  app<-approxSSM(model)
  -sum(out$alpha[,14]*(model$y[,1]-exp(out$theta[,1])))
}



grad(likfn,1)


likfn3<-function(pars){
  model$Z[1,14,]<-pars
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  logLik(app)+ KFAS:::scaling(model$y, model$u, model$distribution, app$theta)-
    sum(dnorm(app$y,app$theta,sqrt(c(t(apply(app$H,3,diag)))),log=T))
}

grad(likfn3,1)

likfn4<-function(pars){  
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  app$Z[1,14,]<-pars  
  out<-KFS(app,smoothing="mean")
  (logLik(app)+ KFAS:::scaling(model$y, model$u, model$distribution, out$mu)-
    sum(dnorm(app$y,out$mu,sqrt(c(t(apply(app$H,3,diag)))),log=T)))
}

grad(likfn4,1)
model2<-model
model2$Z[1,14,]<-1
app2<-approxSSM(model2,tol=1e-15,maxiter=500)
logLik(app2)
x<-sapply(1:37,function(i)
  (app$y[i,1]*out$alpha[i,14]-(out$V[14,,i]+out$alpha[i,14]%*%t(out$alpha[i,]))%*%model$Z[1,,i])/app$H[1,1,i]
) ### sama kuin likfn4!!


######################################
likfn2<-function(pars){
  model$Z[1,14,]<-pars
  logLik(approxSSM(model))
}

likfn2b<-function(pars){
  model$Z[1,14,]<-pars
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  logLik(model,convtol=1e-15,maxiter=500)-KFAS:::scaling(model$y, model$u, model$distribution, app$theta)+
    sum(dnorm(app$y,app$theta,sqrt(c(t(apply(app$H,3,diag)))),log=T))
}
likfn3<-function(pars){
  model$Z[1,14,]<-pars
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  logLik(app)+ KFAS:::scaling(model$y, model$u, model$distribution, app$theta)-
    sum(dnorm(app$y,app$theta,sqrt(c(t(apply(app$H,3,diag)))),log=T))
}

likfnP<-function(pars){
  model$Z[1,14,]<-pars
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  KFAS:::scaling(model$y, model$u, model$distribution, app$theta)
}
likfnG<-function(pars){
  model$Z[1,14,]<-pars
  app<-approxSSM(model,tol=1e-15,maxiter=500)
    sum(dnorm(app$y,app$theta,sqrt(c(t(apply(app$H,3,diag)))),log=T))
}

likfn4<-function(pars){
  app<-approxSSM(model)
  app$Z[1,14,]<-pars
  logLik(app)
}

likfnP2<-function(pars){
  
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  app$Z[1,14,]<-pars
  out<-KFS(app,smoothing="mean")
  KFAS:::scaling(model$y, model$u, model$distribution, out$mu)
}
grad(likfnP,model$Z[1,14,1])
grad(likfnP2,model$Z[1,14,1])

grad(likfn1,model$Z[1,14,1]+1)
grad(likfn2,model$Z[1,14,1])
grad(likfn2,model$Z[1,14,1])
grad(likfn3,model$Z[1,14,1])
grad(likfn4,model$Z[1,14,1])
grad(likfnP,model$Z[1,14,1])
grad(likfnG,model$Z[1,14,1])



likfnP<-function(pars){
  model$Z[1,14,]<-pars
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  sum(dpois(model$y,exp(app$theta),log=TRUE))
}

likfnP2<-function(pars){
  
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  out<-KFS(app)
  app$Z[1,14,]<-pars
  theta<-matrix(0,37,13)
  for(i in 1:37)
    theta[i,]<-app$Z[,,i]%*%out$alpha[i,]
  sum(dpois(model$y,exp(theta),log=TRUE))
}

likfnP3<-function(pars){ 
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  app$Z[1,14,]<-pars
  out<-KFS(app)
  theta<-matrix(0,37,13)
  for(i in 1:37)
    theta[i,]<-app$Z[,,i]%*%out$alpha[i,]
  sum(dpois(model$y,exp(theta),log=TRUE))
}
grad(likfnP,model$Z[1,14,1]+1)
grad(likfnP2,model$Z[1,14,1]+1)
grad(likfnP3,model$Z[1,14,1]+1)

gradf<-function(pars){
  model$Z[1,14,]<-pars
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  out<-KFS(app,smoothing="state")
  gradmat<-0
  k<-14
  #for(j in 1:13)
  j<-1
  for(i in 1:37){
    gradmat<-gradmat+(model$y[i,j]-exp(app$theta[i,j]))*out$alpha[i,k]#+(model$y[i,j]-exp(app$theta[i,j]))*out$alpha[i,k]
  }
  gradmat
}

gradf<-function(pars){
  model$Z[1,14,]<-pars
  app<-approxSSM(model,tol=1e-15,maxiter=500)
  out<-KFS(app,smoothing="state")
  gradmat<-0
  k<-14
  j<-1
  for(i in 1:37){
    gradmat<-gradmat+0.5*out$alpha[i,k]*(1+model$y[i,j]^2*exp(-app$theta[i,j])-exp(app$theta[i,j]))+ #approx, purettu ytilde ja h auki
                     (model$y[i,j]-exp(app$theta[i,j]))*out$alpha[i,k]                              #poisson p(y|theta)
                     (1/app$H[j,j,i])*(app$y[i,j]*out$alpha[i,k]-(out$alpha[i,k]%*%t(out$alpha[i,]))%*%model$Z[j,,i])                                                                            #gaussian g(y|theta) 
                  
  }
  gradmat
}
gradf(model$Z[1,14,1])

gradmat<-0
j<-1
k<-14
for(i in 1:37){
  gradmat<-gradmat+(model$y[i,j]-exp(app$theta[i,j]))*out$alpha[i,k]
}
gradmat


f<-function(pars){
  model$Z[1,14,]<-pars
  app<-approxSSM(model)
  out<-KFS(app)
  x<-0
  for(i in 1:37)
   x<- x + dpois(model$y[i,1],exp(out$mu[i,1]),log=TRUE)
  x
}

f2<-function(pars){
  model$Z[1,14,]<-pars
  out<-KFS(approxSSM(model))
  app<-approxSSM(model)
  x<-0
  for(i in 1:37)
    x<- x + dpois(model$y[i,1],exp(app$theta[i,1]),log=TRUE)
  x
}

grad(likfnP,model$Z[1,14,1])
grad(f,model$Z[1,14,1])
grad(f2,model$Z[1,14,1])
out<-KFS(approxSSM(model))
app<-approxSSM(model)
x<-0
for(i in 1:37)
  x<- x + dpois(model$y[i,1],exp(app$theta[i,1]),log=TRUE)*exp(app$theta[i,1])*out$alpha[i,14]
x
gradmat<-0
j<-1
k<-14
for(i in 1:37){
  gradmat<-gradmat+(model$y[i,j]*out$alpha[i,k]-exp(app$theta[i,j]))
}
gradmat

gradmat<-0
j<-1
k<-14
    for(i in 1:37){
      gradmat<-gradmat+(model$y[i,j]-exp(model$Z[j,,i]%*%out$alpha[i,]))*out$alpha[i,k]
    }
gradmat

x<-sapply(1:37,function(i)
  (app$y[i,1]*out$alpha[i,14]-(out$V[14,,i]+out$alpha[i,14]%*%t(out$alpha[i,]))%*%model$Z[1,,i])/app$H[1,1,i]
) ### sama kuin likfn4!!

sum(x)

likfn2<-function(pars){
  model$Z[1,14,]<-pars
  logLik(approxSSM(model))
}

likfn4<-function(pars){
  app<-approxSSM(model)
  app$Z[1,14,]<-pars
  logLik(app)
}



grad(likfn2,model$Z[1,14,1])/
grad(likfn4,model$Z[1,14,1])



x<-sapply(1:37,function(i)
  (app$y[i,1]*out$alpha[i,14]-(out$V[14,,i]+out$alpha[i,14]%*%t(out$alpha[i,]))%*%model$Z[1,,i])/app$H[1,1,i]*
    out$alpha[i,14]
) ### sama kuin likfn4!!

sum(x)


x<-lapply(1:37,function(i)
  (t(app$y[i,]%*%t(out$alpha[i,]))-(out$V[,,i]+out$alpha[i,]%*%t(out$alpha[i,]))%*%t(model$Z[,,i]))%*%solve(app$H[,,i])
)


(app$y[i,1]*out$alpha[i,14]-(out$V[14,,i]+out$alpha[i,14]%*%t(out$alpha[i,]))%*%model$Z[1,,i])/app$H[1,1,i]

x<-sapply(1:37,function(i)
(app$theta[i,1]*out$alpha[i,14]-(out$V[14,,i]+out$alpha[i,14]%*%t(out$alpha[i,]))%*%model$Z[1,,i])/app$H[1,1,i]
)


x<-sapply(1:37,function(i)
  (app$theta[i,1]*out$alpha[i,]-(out$V[,,i]+out$alpha[i,]%*%t(out$alpha[i,]))%*%model$Z[1,,i])/app$H[1,1,i]
)

sum(x)
gradmat1<-gradmat2<-gradmat3<-matrix(0,13,2)

for(j in 1:13)
  for(k in 14:15)
    for(i in 1:37){
      gradmat1[j,k-13]<-gradmat1[j,k-13]+(app$y[i,j]*out$alpha[i,k]-(out$V[k,,i]+out$alpha[i,k]%*%t(out$alpha[i,]))%*%model$Z[j,,i])/app$H[j,j,i]
      gradmat2[j,k-13]<-gradmat2[j,k-13]+(model$y[i,j]-exp(model$Z[j,,i]%*%out$alpha[i,]))*out$alpha[i,k]
      gradmat3[j,k-13]<-
        gradmat3[j,k-13]+(1/(2*app$H[j,j,i]))*(-2*app$y[i,j]*out$alpha[i,k]+2*(out$alpha[i,k]%*%t(out$alpha[i,]))%*%(model$Z[j,,i]))    
    }
gradmat1[1,1]
gradmat2[1,1]
gradmat3[1,1]
library(numDeriv)


out0<-logLik(fit_glmm$model,maxit=0,smoothing=c("state","signal","mean"),simp=F)
outa0<-KFS(approxSSM(fit_glmm$model,maxit=0),smoothing=c("state","signal","mean"),simp=F)
out0
outa0$log

out1<-logLik(fit_glmm$model,maxit=1,smoothing=c("state","signal","mean"),simp=F)
outa1<-KFS(approxSSM(fit_glmm$model,maxit=1),smoothing=c("state","signal","mean"),simp=F)
out1
outa1$log
unname(outa0$alpha[1,])

[1] -36.513980 -24.545513  10.868912  -6.640139  -6.094058 -15.179102 -25.905942  11.640785  72.814924 -28.219073
[11] -30.123995 105.197566 -44.406528  36.424904  43.994216  13.805380 -13.083712  30.577154  31.069692  12.954167
[21]  -7.898943  23.036976  48.287925 103.192641  21.590591

fit_glmm<-glmm(response.var="abundance",grouping.var="species",fixed=~1,maxiter=100,convtol=1e-10,
               latent.factors=2,init.factors=initmat_glmm[lower.tri(initmat_glmm,TRUE)],
               data=borneodata,distribution="poisson",print_level=1,estimate=F)
likfn<-function(pars,model){
  model$Z[,14:15,]<-c(pars[1:13],0,pars[14:25])
  -logLik(model)
}
likfn2<-function(pars,model){
  model$Z[,14:15,]<-c(pars[1:13],0,pars[14:25])
  apx<-approxSSM(model)
  -(logLik(model)- KFAS:::scaling(model$y, model$u, model$distribution, apx$theta)+sum(dnorm(apx$y,apx$theta,sqrt(c(t(apply(apx$H,3,diag)))),log=T)))
}

likfn3<-function(pars,model){
  model$Z[,14:15,]<-c(pars[1:13],0,pars[14:25])
  apx<-approxSSM(model)
  -logLik(apx)
}


pars<-initmat_glmm[lower.tri(initmat_glmm,TRUE)]
grad(likfn,pars,model=fit_glmm$model)
gg2<-grad(likfn3,pars,model=fit_glmm$model)
grad(likfn2,pars,model=fit_glmm$model)
model<-fit_glmm$model

model$Z[,14:15,]<-c(pars[1:13],0,pars[14:25])
app<-approxSSM(model) 
out<-KFS(app,smoothing="state")
gradmat1<-gradmat2<-gradmat3<-matrix(0,13,2)

for(j in 1:13)
  for(k in 14:15)
   for(i in 1:37){
  gradmat1[j,k-13]<-gradmat1[j,k-13]+(app$y[i,j]*out$alpha[i,k]-(out$V[k,,i]+out$alpha[i,k]%*%t(out$alpha[i,]))%*%model$Z[j,,i])/app$H[j,j,i]
  gradmat2[j,k-13]<-gradmat2[j,k-13]+(model$y[i,j]-exp(model$Z[j,,i]%*%out$alpha[i,]))*out$alpha[i,k]
  gradmat3[j,k-13]<-
    gradmat3[j,k-13]+(1/(2*app$H[j,j,i]))*(-2*app$y[i,j]*out$alpha[i,k]+2*(out$alpha[i,k]%*%t(out$alpha[i,]))%*%(model$Z[j,,i]))    
}
-c(gradmat1)
-c(gradmat2)
-c(gradmat3)


likfnG<-function(pars,model){
  model$Z[,14:15,]<-
  apx<-approxSSM(model)
  -dnorm(apx$y[1,1],apx$theta[1,1],sqrt(apx$H[1,1,1]),log=T)
}
likfnG<-function(pars,model){
  model$Z[1,14,1]<-pars
  apx<-approxSSM(model)
  -dnorm(apx$y[1,1],apx$theta[1,1],sqrt(apx$H[1,1,1]),log=T)
}

pars<-initmat_glmm[lower.tri(initmat_glmm,TRUE)]
model$Z[,14:15,]<-c(pars[1:13],0,pars[14:25])
pars<-2
grad(likfnG,pars,model=model)
model$Z[1,14,1]<-pars
apx<-approxSSM(model)
out<-KFS(apx,smooth="state")
(-(1/apx$H[1,1,1])*(-apx$y[1,1]*out$alp[1,14]+(out$alp[1,14]%*%t(out$alp[1,])%*%model$Z[1,,1])))/dnorm(apx$y[1,1],apx$theta[1,1],sqrt(apx$H[1,1,1]))


model$Z[,14:15,]<-c(pars[1:13],0,pars[14:25])
model$Z[,14,1]<-pars
app<-approxSSM(model) 
out<-KFS(app,smoothing="state")
j<-1
k<-14
i<-1
1/(app$H[j,j,i])*(-app$y[i,j]*out$alpha[i,k]+(out$alpha[i,k]%*%t(out$alpha[i,])%*%model$Z[j,,i]))
1/(app$H[j,j,i])*(-app$y[i,j]*out$alpha[i,]+(out$alpha[i,]%*%t(out$alpha[i,]))%*%model$Z[j,,i])

c(solve(app$H[,,i])*(-app$y[i,j]*out$alpha[i,]+((out$alpha[i,]%*%t(out$alpha[i,])))%*%model$Z[,,i])

 likfnG<-function(pars,model){
  model$Z[,14:15,]<-c(pars[1:13],0,pars[14:25])
  apx<-approxSSM(model)
  -sum(dnorm(apx$y,apx$theta,sqrt(c(t(apply(apx$H,3,diag)))),log=T))
}

likfnG<-function(pars){
  x<-c(1,pars)
  y<-5
  a<-c(2,3)
  h<-1
  dnorm(y,x%*%a,sqrt(h),log=T)
}

gradf<-function(pars){
  x<-c(1,pars)
  y<-5
  a<-c(2,3)
  h<-1
-(1/h)*(-y*a[2]+(a[2]%*%t(a)%*%x))
}
gradf(5)
grad(likfnG,5)

likfnG<-function(pars){
  x<-c(pars,pars)
  y<-5
  a<-c(2,3)
  h<-1
  dnorm(y,x%*%a,sqrt(h),log=T)
}

gradf<-function(pars){
  x<-c(pars,pars)
  y<-5
  a<-c(2,3)
  h<-1
  -(1/h)*(-y*a[1]+(a[1]%*%t(a)%*%x))-(1/h)*(-y*a[2]+(a[2]%*%t(a)%*%x))
}
gradf(5)
grad(likfnG,5)


tmpm = coefVars.slice(t).row(k);
tmpm2 = coefs.col(t).t();
grad(j,k) += (ytilde(t,j)*coefs(k,t)-
                arma::as_scalar((tmpm.cols(zind.col(j)) + coefs(k,t)*tmpm2.cols(zind.col(j)))*Z.slice(t).col(j)))/H(t,j); 
likfn(initmat_glmm[lower.tri(initmat_glmm,TRUE)],model=fit_glmm$model)

all.equal(unname(fit_glmm$fixed$coef),fit_glmmf1$coefs[3:15,37])
all.equal(unname(fit_glmm$fixed$coef),fit_glmmf1$coefs[3:15,1])
all.equal(t(unname(fit_glmm$factors$u)),fit_glmmf1$coefs[1:2,])
out<-KFS(fit_glmm$model,smoothing=c("state","signal","mean"))

all.equal(t(unclass(unname(out$alp)[,1:13])),fit_glmmf1$coefs[3:15,],check.attributes=FALSE)
unname(out$alp)[,1]
all.equal(unname(t(outa0$alp)),fit_glmmf1$coefs)
unname(outa0$r[1,])
unname(outa0$alp[1,])
-643.041
-1252.48
-0.9510  -0.7267   2.2972   2.0306   1.5586   1.1562   1.3007   1.6199   1.9270   1.3222   2.6679   2.1350   2.4867   2.5127   1.8053

Iteration 1, Log-likelihood: -1513.66

fit_glmmf2<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
                  convtol=1e-8,maxiter=25,maxiter2=0,init.factor=initmat_glmmf,
                  data=borneodata,control=list(trace=1,REPORT=1,maxit=2),method="BFGS",
                  estimate=T,trace=1,gradient=FALSE)
fit_glmmf1$logLik


fit_glmmf$coefs
fit_glmm<-suppressWarnings(glmm(response.var="abundance",grouping.var="species",fixed=~1,maxiter=2,
          latent.factors=2,init.factors=initmat_glmm[lower.tri(initmat_glmm,TRUE)],data=borneodata,distribution="poisson",print_level=1,estimate=FALSE))
fit_glmm$logLik
#-1185.676
#-1392.838
#-1374.671

suppressWarnings(deviance(KFS(glmm(response.var="abundance",grouping.var="species",fixed=~Logging+Slope,maxiter=0,
                                latent.factors=2,init.factors=initmat_glmm[lower.tri(initmat_glmm,TRUE)],
                                data=borneodata,distribution="poisson",print_level=1,estimate=FALSE)$model,maxiter=2)))
#1187.882
#979.3993
#968.1192
#-1277.913



factors<-fit_glmmf$coefs[1:2,]
plot(x=factors[1,],y=factors[2,],pch=NA)
text(x=factors[1,],y=factors[2,], labels=as.character(1:37))

dat = read.csv("born.csv")
Y = dat[,-(1:4)]
X = dat[,1:4] #remove id column
plot(x=factors[1,],y=factors[2,],type="n",main="LV + main effects",xlab="LV1",ylab="LV2")
points(x=factors[1,],y=factors[2,],col=as.numeric(X[,2]),pch=as.numeric(X[,3]))
legend("topleft", c("P","L89","L93","L","M","U"),col=c(3,1,2,1,1,1),pch=c(15,15,15,1,2,3))

#######
app<-approxSSM(fit_glmm$model,maxit=0)
logLik(fit_glmm$model,maxit=0)-logLik(app)
#[1] -609.4362
sum(dpois(c(fit_glmm$model$y),c(exp(app$theta)),log=TRUE)-dnorm(c(app$y),c(app$theta),sqrt(c(t(apply(app$H,3,diag)))),log=TRUE))
#[1] -609.4362

logLik(app)+sum(dpois(c(fit_glmm$model$y),c(exp(app$theta)),log=TRUE)-dnorm(c(app$y),c(app$theta),sqrt(c(t(apply(app$H,3,diag)))),log=TRUE))
-1185.676

dpois(fit_glmm$model$y[1],exp(app$theta[1]),log=TRUE)-dnorm(app$y[1],app$theta[1],sqrt(app$H[1]),log=TRUE)
##-424.9182

app<-approxSSM(fit_glmm$model,maxit=1)


out<-KFS(approxSSM(fit_glmm$model,maxit=0),filtering=c("state","mean"),smoothing=c("state","mean"),simp=F)
out$F[,1]
unname(out$a[2,])
fit_glmmf$coefs[,1]

unname(out$a[38,])
fit_glmmf$coefs[,37]

all.equal(fit_glmmf$lin,unclass(out$mu),check.attributes=FALSE)
logLik(fit_glmm$model,maxit=0)
fit_glmmf$logLik #-1462.26


unname(fit_glmm$fixed$coef)
fit_glmmf$coefs[,37]
out<-KFS(approxSSM(fit_glmm$model,maxiter=1),filtering="state")
unname(out$a[2,])
fit_glmmf$coefs[,1]
unname(out$a[38,])
fit_glmmf$coefs[,37]


fit_glmmf$lin[1:10]
app$theta[1:10]



fit_glmmf<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
           convtol=1e-8,maxiter=1,maxiter2=0,init.theta=fitinit_glmmf$linear.predictor,
           data=borneodata,control=list(trace=1,REPORT=1,maxit=2),method="BFGS",init.factor=initmat_glmmf,estimate=F,trace=2)
fit_glmmf$logLik

initmat<-t(matrix(fitinit$fixed$coef,ncol=13)[6:7,])

init.factors<-initmat[lower.tri(initmat,TRUE)]

a<-proc.time()


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
