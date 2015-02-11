
set.seed(1)
randu1<-rnorm(37)
randu2<-rnorm(37)

dataf<-cbind(borneodata,randu1,randu2)
save(dataf,file="dataf.rda")

model<-  buildGLMMF(group="species",response="abundance",distinct=~randu1+randu2,data=dataf,distribution="poisson")

f<-function(maxiter){
theta <- KFAS:::init_theta(model$y, model$u, model$distribution) 
dist<-pmatch(x = model$distribution, 
             table = c("gaussian", "poisson", "binomial", 
                       "gamma", "negative binomial"), duplicates.ok = TRUE)

res<-GLMMF:::expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                             dist, model$tol, maxiter, 1e-10, theta, model$Zind, model$nfactors)


res$logLik+KFAS:::scaling(model$y,model$u,model$distribution,res$theta)
}
f(100) ##-9712.146
microbenchmark(f(100),times=5) #5.34463

model<-  buildGLMMF(group="species",response="abundance",distinct=~Logging+Slope+randu1+randu2,data=dataf,distribution="poisson")
length(model$a1) #959 tilaa
dim(model$Z) #7 137  37

f(maxiter=0) 
microbenchmark(f(0),times=1) #5.556405
microbenchmark(f(5),times=1) #28.0284

microbenchmark(expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                               dist, model$tol,1,1e-10,theta),times=2)
