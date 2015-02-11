likfn<-function(pars,model,estimate=TRUE,env=parent.frame()){     
  
  if(nf>0){
    Z<-matrix(0,nrow=p,ncol=nf)
    Z[lower.tri(Z,TRUE)]<-pars[1:length(pars)]
    model$Z[1:nf,1:p,]<-t(Z)
  }
  if(estimate){        
    if(itertheta)
      theta <-get("prevtheta",envir=env)
    
    res<-GLMMF:::expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                         dist, model$tol, maxiter, maxiter2, convtol, theta, model$Zind, model$nfactors,trace,0)
    res$theta[res$theta>5]<-5
    res$theta[res$theta< -5]<- -5
    assign("prevtheta",res$theta,envir=env)
    if(res$conv<0){
      warning("Approximating algorithm did not converge.")
      lik<- .Machine$double.xmax                  
    } else lik <- -res$logLik
    lik
  }else model
}

env<-new.env()

nf<-2
p<-ncol(model$Z)
itertheta<-TRUE
maxiter<-100
maxiter2<-0
convtol<-1e-10
dist<-2
trace<-2
library(numDeriv)
pars<-fitBFGS$fit$par
itertheta<-TRUE
env$prevtheta<-fitBFGS$linear.predictor
numgrad1<-grad(f=likfn,x=pars,model=model,env=env)

env$prevtheta<-fitBFGS$linear.predictor
numgrad2<-grad(f=likfn,x=pars,model=model,env=env)

numgrad3<-grad(f=likfn,x=pars,model=model,env=env)
all.equal(numgrad1,numgrad2)
all.equal(numgrad1,numgrad3)


grfn<-function(pars,model,estimate=TRUE,env=parent.frame()){     
  
  if(nf>0){
    Z<-matrix(0,nrow=p,ncol=nf)
    Z[lower.tri(Z,TRUE)]<-pars[1:length(pars)]
    model$Z[1:nf,1:p,]<-t(Z)
  }
  if(estimate){        
    if(itertheta)
      theta <-get("prevtheta",envir=env)
    
    res<-GLMMF:::expfLogLikNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                                 dist, model$tol, maxiter, maxiter2, convtol, theta, model$Zind, model$nfactors,trace,1)
    res$theta[res$theta>5]<-5
    res$theta[res$theta< -5]<- -5
    assign("prevtheta",res$theta,envir=env)
    if(res$conv<0){
      warning("Approximating algorithm did not converge.")
      lik<- 0               
    } else lik <- -res$gradient
    lik
  }else model
}

grfn(pars,model=model,env=env)

theta<-fitBFGS$linear.predictor
out<-GLMMF:::filterNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, dist, 
                 model$tol, maxiter, maxiter2, convtol, theta, model$Zind, model$nfactors,trace)

gradmat<-matrix(0,13,2)
H <- exp(-out$linear.predictor)
ytilde<-model$y*H + out$linear.predictor -1


for(t in 1:37)
  for(j in 1:13)
  gradmat[j,] <- (model$y[t,j]*out$coefs[1:2,t]-(out$coefVars[1:2,1:2]+out$coefs[1:2,t]%*%t(out$coefs[1:2,t]))*model$Z[1])/H[t,j]