load("../StateSpaceGLMM/borneodata.rda")

dat = read.csv("../StateSpaceGLMM/born.csv")
Y = dat[,-(1:4)]
X = dat[,2:4] #remove id column
#remove rare species:
Y = as.matrix(Y[,apply(Y>0,2,sum)>30])

set.seed(1)
randu1<-rnorm(37)
randu2<-rnorm(37)

save(Y,dat,X,file="D:/dat.rda")

fitinit<-glmmf(group="species",response="abundance",distinct=~Logging+Slope+randu1+randu2,
               nfactors=0,distribution="poisson",
               data=borneodata,control=list(trace=1,REPORT=1),method="BFGS")

fitinit$log #-3307.467

initmat<-t(matrix(fitinit$coefs[,28],nrow=3)[2:3,])



fit<-glmmf(group="species",response="abundance",distinct=~1,nfactors=2,distribution="poisson",
           data=spiderdata,control=list(trace=1,REPORT=1,maxit=1000),method="BFGS",init.factor=initmat,estimate=TRUE)

factors<-fit$coefs[1:2,]
plot(x=factors[1,],y=factors[2,],pch=NA)
text(x=factors[1,],y=factors[2,], labels=as.character(1:28))



model<-  buildSSGLMM(Y ~ randu1+randu2,distribution="poisson")

f<-function(maxiter=100,convtol=1e-10){
  dist<-pmatch(x = model$distribution, 
               table = c("gaussian", "poisson", "binomial", 
                         "gamma", "negative binomial"), duplicates.ok = TRUE)
  theta<-KFAS:::init_theta(model$y, model$u, model$distribution)
  
  out<-expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                       dist, model$tol,maxiter,convtol,theta)
  out$logLik+KFAS:::scaling(model$y,model$u,model$distribution,out$theta)
}
f() ##-9712.146
microbenchmark(f(),times=1) #13.50376

model<-  buildSSGLMM(Y ~ Logging+Slope+randu1+randu2,data=X,distribution="poisson")
dim(model$a1) #959 tilaa
dim(model$Z) #7 137  37

f(maxiter=0) 
microbenchmark(f(),times=1)


microbenchmark(expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                     dist, model$tol,1,1e-10,theta),times=2)

