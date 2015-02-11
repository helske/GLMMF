library(KFAS)
library(RcppSSGLMM)

set.seed(1)
counts<-matrix(sample(1:30,replace=TRUE,size=1000),ncol=10)
counts[sample(1:length(counts),size=100)]<-NA
outcome <- gl(5,1,100)
treatment <- gl(5,20)
model<-buildSSGLMM(counts ~ outcome + treatment,distribution="gaussian")

modelKFAS<-SSModel(counts ~ outcome + treatment,distribution="gaussian")

logLik(modelKFAS) #-4396.579

f<-function(maxiter){
  dist<-pmatch(x = model$distribution, 
               table = c("gaussian", "poisson", "binomial", 
                         "gamma", "negative binomial"), duplicates.ok = TRUE)
  theta<-KFAS:::init_theta(model$y, model$u, model$distribution)
  
  out<-expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                       dist, model$tol,maxiter,1e-15,theta)
  out$logLik+KFAS:::scaling(model$y,model$u,model$distribution,out$theta)
}

diag(modelKFAS$P1)<-1
diag(modelKFAS$P1inf)<-0
diag(model$P1)<-1
diag(model$P1inf)<-0
f(0)
logLik(modelKFAS,maxiter=0)

system.time(f())
f() #-4396.579
microbenchmark(logLik(modelKFAS),f(),times=3)

####
counts <- c(18,17,15,20,10,20,25,13,12)
outcome <- gl(3,1,9)
treatment <- gl(3,3)
print(d.AD <- data.frame(treatment, outcome, counts))
glm.D93 <- glm(counts ~ outcome + treatment, family = poisson())
summary(glm.D93)

model<-buildSSGLMM(counts ~ outcome + treatment,distribution="poisson")

dist<-pmatch(x = model$distribution, 
             table = c("gaussian", "poisson", "binomial", 
                       "gamma", "negative binomial"), duplicates.ok = TRUE)

f<-function(){
theta<-KFAS:::init_theta(model$y, model$u, model$distribution)
expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                dist, model$tol,100,1e-10,theta)
}

modelKFAS<-SSModel(counts ~ outcome + treatment,distribution="poisson")
app<-approxSSM(modelKFAS,maxiter=100)
app$y
app$H
###
counts<-matrix(sample(1:30,replace=TRUE,size=10000),ncol=20)
counts[sample(1:length(counts),size=1000)]<-NA
outcome <- gl(5,1,500)
treatment <- gl(5,100)
model<-buildSSGLMM(counts ~ outcome + treatment,distribution="poisson")

modelKFAS<-SSModel(counts ~ outcome + treatment,distribution="poisson")

logLik(modelKFAS)
f<-function(){
theta<-KFAS:::init_theta(model$y, model$u, model$distribution)
dist<-pmatch(x = model$distribution, 
             table = c("gaussian", "poisson", "binomial", 
                       "gamma", "negative binomial"), duplicates.ok = TRUE)
lik<-expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                     dist, model$tol,100,1e-10,theta)
lik+KFAS:::scaling(model$y,model$u,model$distribution,theta)
}
microbenchmark(logLik(modelKFAS),f(),times=3)
#Unit: seconds
#expr      min       lq     mean   median       uq      max neval cld
#logLik(modelKFAS) 60.46974 60.51758 60.72329 60.56542 60.85007 61.13471     3   b
#f()               12.43671 12.44755 12.46990 12.45839 12.48650 12.51461     3  a 

dist<-pmatch(x = model$distribution, 
             table = c("gaussian", "poisson", "binomial", 
                       "gamma", "negative binomial"), duplicates.ok = TRUE)

f<-function(){
  theta<-KFAS:::init_theta(model$y, model$u, model$distribution)
  lik<-expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                  dist, model$tol,100,1e-10,theta)
  lik+KFAS:::scaling(model$y,model$u,model$distribution,theta)
}

modelKFAS<-SSModel(counts ~ outcome + treatment,distribution="poisson")
app<-approxSSM(modelKFAS,maxiter=100)

all.equal(logLik(modelKFAS),f())
microbenchmark(logLik(modelKFAS),f())
logLik(app)
out<-KFS(app)
(sum(dnorm(c(app$y),mean=c(app$theta),sd=sqrt(c(t(apply(app$H,3,diag)))),log=TRUE))-0.5*10*log(2*pi))

logLik(modelKFAS)
f()-sum(dnorm(c(app$y),mean=c(app$theta),sd=sqrt(c(t(apply(app$H,3,diag)))),log=TRUE))
#####
ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
weight <- c(ctl, trt)
lm.D9 <- lm(weight ~ group)

model<-SSModel(weight ~ group)
logLik(model)
o<-KFS(model,smo="mean")
(logLik(model)-sum(dnorm(c(model$y),mean=c(o$mu),sd=1,log=TRUE)))/(0.5*log(2*pi))
+2*log(2*pi)


out<-kfilter(model$y, model$Z, model$H, model$a1, model$P1, model$P1inf, TRUE,model$tol)
model$H[]<-sum((model$y-out$fitted)^2)/18
out2<-kfilter(model$y, model$Z, model$H, model$a1, model$P1, model$P1inf, TRUE,model$tol)
out2$logLik
logLik(lm.D9,REML=TRUE)

model2<-model
model2$H[]<-NA
fit<-fitSSM(model=model2,inits=0,method="BFGS")

microbenchmark(glmmSSLogLikGaussian(model$y, model$Z, model$H, model$a1, model$P1, model$P1inf, 1e-8),logLik(model))

microbenchmark(logLik(model),logLik2(model)) #loglik2:ssa glmSSlogLik sijoitettu logLikiin.
#Unit: microseconds
#expr     min      lq     mean   median     uq     max neval cld
#logLik(model) 431.153 464.089 491.2710 469.2210 491.25 977.793   100   b
#logLik2(model) 417.038 442.488 459.7471 448.9045 465.80 853.752   100  a 

y<-matrix(rnorm(3000),100,30)
x1<-rnorm(100)
x2<-rnorm(100)
x3<-rnorm(100)
x4<-rnorm(100)
x5<-rnorm(100)
x6<-rnorm(100)
m<-SSModel(y~x1+x2+x3+x4+x5+x6)

logLik(m)+0.5*3000*log(2*pi)
logLik2(m)
microbenchmark(logLik(m),logLik2(m),times=5)
logLik(model)+0.5*18*log(2*pi)-l

model2<-model
model2$P1inf[]<-0
diag(model2$P1)<-1e7
logLik(model2)+0.5*20*log(2*pi)
glmmSSLogLikGaussian(model2$y, model2$Z, model2$H, model2$a1, model2$P1, model2$P1inf, model2$tol)
