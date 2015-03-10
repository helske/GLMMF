fitGLMMF2<-function(model,estimate.dispersion,correlating.effects,maxiter,maxiter2, convtol,
                   init.random.cov,init.dispersion,init.factor,init.theta=NULL,
                   common.dispersion,estimate,trace,
                   gradient=TRUE,maxapp=10,epsilon=1e-5,...){
  
  
  #for more compact likfn  
  nr<-model$nrandom
  nf<-model$nfactors
  nc<-model$ncommon
  nd<-model$ndistinct
  n<-nrow(model$y)
  p<-ncol(model$y)
  m<-nrow(model$P1)
  if(!is.null(init.theta) && !identical(dim(init.theta),c(n,p))){    
    stop("init.theta is not of valid form.")       
  } 
  
  if(missing(estimate.dispersion) || estimate.dispersion){
    if(model$distribution%in%c("gaussian","negative binomial","gamma")){
      ndisp<-ifelse(common.dispersion,1,p)
    } else ndisp<-0
  } else ndisp<-0
  
  if(nr<2) correlating.effects<-FALSE
  
  if(nf>0){
    if(missing(init.factor)){
      init.factors<-rep(1,p*nf-sum(upper.tri(diag(nf))))
    } else {
      if(!identical(dim(init.factor),c(p,as.integer(nf)))){
        stop("Dimensions of the initial latent factor matrix is not equal to p*nfactors.")
      }
      init.factors<-init.factor[lower.tri(init.factor,TRUE)]
    }  
  } else init.factors<-NULL
  
  if(ndisp>0){
    if(missing(init.dispersion)){
      init.dispersion<-rep(1,ndisp)
    } else {
      if(length(init.dispersion)==1){
        init.dispersion<-rep(log(init.dispersion),ndisp)
      } else {
        if(length(init.dispersion)!=ndisp) 
          stop("Number of initial values for distinct dispersion parameters is not equal to groups")
        init.dispersion<-log(init.dispersion)
      }
    }  
  }
  
  if(nr>0){
    if(missing(init.random.cov)){
      if(correlating.effects){ 
        init.random <- diag(nr)[upper.tri(diag(nr),TRUE)] 
      }else{ 
        init.random <- rep(1, nr)
      }
    }else {
      if(!identical(dim(init.random.cov),c(nr,nr)))
        stop("Dimensions of the initial covariance matrix of random effects does not match with number of random effects.")
      
      if(correlating.effects){
        init.random.cov<-VarCorr(fm1)$Subject
        tmp<-chol(init.random.cov)
        init.random<-c(2*log(diag(tmp)),tmp[upper.tri(tmp)])  
      } else init.random<-sqrt(init.random.cov[1 + 0:(nr - 1) * (nr + 1)])
    }
  } else init.random<-NULL

  if(is.null(init.theta)){
    theta <- init_theta(model$y, model$u, model$distribution) 
  } else {
    theta <- init.theta
  }
  
 
  

  if(correlating.effects){
    likfn<-function(pars,model,estimate=TRUE){
#       if(ndisp>0){        
#         model$u[]<-matrix(exp(pars[1:ndisp]),n,p,byrow=TRUE)
#       }
      if(nr>0){
        if(any(exp(0.5*pars[ndisp+1:nr])>1e5)) #mixing with fixed effects
          return(.Machine$double.xmax)
        P1<-diag(exp(0.5*pars[ndisp+1:nr]))
        P1[upper.tri(P1)]<-pars[ndisp+nr+1:length(init.random[-(1:nr)])]
        P1<-crossprod(P1)
        model$P1[(nf+nc+nd+1):m,(nf+nc+nd+1):m]<-
          as.matrix(.bdiag(replicate(p,P1,simplify=FALSE)))
      }
      if(nf>0){
        Z<-matrix(0,nrow=p,ncol=nf)
        Z[lower.tri(Z,TRUE)]<-pars[(ndisp+(nr>0)*length(init.random)+1):length(pars)]
        model$Z[1:nf,1:p,]<-t(Z)
      }
      if(estimate){
#         if(itertheta)
#           theta <-get("prevtheta",pos=env)
        
        res<-likelihood(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                             1, model$tol, maxiter, maxiter2, convtol, theta, model$Zind, 
                             model$nfactors,trace,1)
#         res$theta[res$theta>5]<-5
#         res$theta[res$theta< -5]<- -5
#         assign("prevtheta",res$theta,pos=env)
        if(res$conv<0){
          warning("Approximating algorithm did not converge.")
          lik<- -.Machine$double.xmax
          if(gradient)
          attr(lik,"gradient")<-0
        } else {
          lik <- res$logLik
          if(gradient)
          attr(lik,"gradient")<- -res$gradient[lower.tri(res$gradient,TRUE)]
        }
        -lik
      }else model
    }
    inits<-c(if(ndisp>0) init.dispersion else NULL, if(nr>0) init.random else NULL,
             if(nf>0) init.factors else NULL)
    
  } else {
    likfn<-function(pars,model,estimate=TRUE){      
#       if(ndisp>0){        
#         model$u[]<-matrix(exp(pars[1:ndisp]),n,p,byrow=TRUE)
#       }
      if(nr>0){
        if(any(exp(pars[ndisp+1:nr])>1e5)) #mixing with fixed effects
          return(.Machine$double.xmax)
        P1<-diag(exp(pars[ndisp+1:nr]))
        
        model$P1[(nf+nc+nd+1):m,(nf+nc+nd+1):m]<-
          as.matrix(.bdiag(replicate(p,P1,simplify=FALSE)))
      }
      if(nf>0){
        Z<-matrix(0,nrow=p,ncol=nf)
        Z[lower.tri(Z,TRUE)]<-pars[(ndisp+nr+1):length(pars)]
        model$Z[1:nf,1:p,]<-t(Z)
      }
      if(estimate){        
#           if(itertheta)
#             theta <-get("prevtheta",envir=env)
#          
        res<-likelihood(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                             1, model$tol, maxiter, maxiter2, convtol, theta, model$Zind, model$nfactors,trace,1)
#         res$theta[res$theta>5]<-5
#         res$theta[res$theta< -5]<- -5
#         assign("prevtheta",res$theta,envir=env)
        if(res$conv<0){
          warning("Approximating algorithm did not converge.")
          lik<- -.Machine$double.xmax
          if(gradient)
           attr(lik,"gradient")<-0
        } else {
          lik <- res$logLik
          if(gradient)
            attr(lik,"gradient")<- -res$gradient[lower.tri(res$gradient,TRUE)]
        }
        -lik
      }else model
    }
    inits<-c(if(ndisp>0) init.dispersion else NULL, if(nr>0) init.random else NULL,
             if(nf>0) init.factors else NULL)  
    
  }
#   env<-new.env()
#   env$prevtheta<-theta
  if(length(inits)>0 && estimate){
    
    amodel<-model         
    dist<-pmatch(x = model$distribution, 
                 table = c("gaussian", "poisson", "binomial", 
                           "gamma", "negative binomial"))    
    true.model<-likelihood(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                                dist, model$tol, maxiter, maxiter2, convtol, theta, model$Zind, model$nfactors,trace,0)     
    theta<-true.model$theta
    llold<-true.model$logLik
    ii<-0
    ll<-10*llold
    while(ii<maxapp && abs(ll - llold)/(0.1 + abs(ll)) > epsilon){
      llold<-ll
      approx.model<-approxG(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                            dist, model$tol, maxiter, maxiter2, convtol, theta, model$Zind, model$nfactors,trace)
      theta<-approx.model$theta
      amodel$u[]<-approx.model$H
      amodel$y[]<-approx.model$y      
      fit<-nlm(p=inits,f=likfn,model=amodel,estimate=TRUE,...)   
      inits<-fit$e
      amodel<-likfn(inits,amodel,FALSE)
      model<-likfn(inits,model,FALSE)      
      true.model<-likelihood(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, 
                      dist, model$tol, maxiter, maxiter2, convtol, theta, model$Zind, model$nfactors,trace,0)     
      theta<-true.model$theta
      ll<-true.model$logLik
      ii <- ii+1
      if(trace>0)
        print(paste("ApproxMaximization iter",ii,"log-likelihood",ll))
    }    
   
    
  } else{
    model<-likfn(inits,model,FALSE)
    fit<-"Initial model returned."
  }
  list(model=model,fit=fit)
}