#' Build Generalized Linear Mixed-Effects Model with Latent Factors
#'
#' Function \code{buildGLMMF} builds GLMM with latent factors used by \code{glmmf}. This function is exported for debugging purposes.
#'
#' @param group Name of the grouping variable in data. 
#' Only one grouping variable is allowed and the group sizes must be equal. In case of unequal group sizes, 
#' patch missing rows with NA's.
#' @param response Name of the response variable in data.
#' @param common.fixed formula for common fixed effects. LHS of the formula is ignored if present.
#' @param distinct.fixed Formula for distinct fixed effects i.e. each group has separate regression 
#' coefficient. This formula cannot contain variables which are already present in 
#' \code{common.fixed} or \code{random} formulas, as in that case the model would not be identifiable. 
#' LHS of the formula is ignored if present.
#' @param random Formula for random effects. LHS of the formula is ignored if present.
#' @param data Data frame containing the variables in the model. Must contain all variables used formulas and 
#' variables defined in \code{group} and \code{response}
#' @param distribution Character. Possible values are "gaussian", "poisson", 
#' "binomial", "negative binomial" and "gamma". Default is "gaussian".
buildGLMMF<-
  function(group,response,common.fixed,distinct.fixed,random,nfactors=0,data,
           distribution=c("gaussian", "poisson", "binomial", "gamma", "negative binomial"),
           u,random.cov,factors,tol=.Machine$double.eps^0.5){
    
    distribution<-match.arg(distribution)
    
    
    splitdata<-split(data,factor(data[,group]))
    if(length(unique(sapply(splitdata,nrow))) != 1)
      stop("Need balanced group sizes. Patch data frame with NAs.")
   
    model<-NULL
    
    model$call <- match.call(expand.dots = FALSE)
    
    
    if(!missing(common.fixed)){
      common_mm<-sapply(1:length(splitdata), function(i) model.matrix(common.fixed,data=splitdata[[i]]),simplify="array")
      ncommon<-dim(common_mm)[2]
      model$terms.common <- terms(common.fixed, data = data)
    } else {
      ncommon<-0
      common_mm<-NULL
    }
    if(!missing(random)){
      random_mm<-sapply(1:length(splitdata), function(i) model.matrix(random,data=splitdata[[i]]),simplify="array") 
      nrandom<-dim(random_mm)[2]
      model$terms.random <- terms(random, data = data)
    } else {
      nrandom <-0
      random_mm <- NULL
    }
    if(!missing(distinct.fixed)){
      distinct_mm<-sapply(1:length(splitdata), function(i) model.matrix(distinct.fixed,data=splitdata[[i]]),simplify="array")
      ndistinct<-dim(distinct_mm)[2]
      if(!missing(common.fixed) && attr(terms(distinct.fixed),"term.labels")%in% attr(terms(common.fixed),"term.labels"))
        stop("Model is not identifiable, variables in distinct.fixed also present in  common.fixed formula.")
      if(!missing(random) && attr(terms(distinct.fixed),"term.labels")%in% attr(terms(random),"term.labels"))
        stop("Model is not identifiable, variables in distinct.fixed also present in random formula.")
      model$terms.distinct <- terms(distinct.fixed, data = data)
    } else {
      ndistinct<-0
      distinct_mm<-NULL
    }
    
    if(nfactors>0){
      if(missing(factors)){
        factors<-matrix(NA,length(splitdata),nfactors)
        factors[upper.tri(factors)]<-0  
      } else {
        if(!identical(dim(factors),c(length(splitdata),as.integer(nfactors)))){
          stop("Dimensions of the latent factor matrix is not equal to p*nfactors.")
        }
      }
    } else factors<-NULL
    
    y<-sapply(1:length(splitdata), function(i) splitdata[[i]][,response],simplify="matrix")
    storage.mode(y)<-"double"
    colnames(y)<-names(splitdata)
    
    p<-ncol(y)
    n<-nrow(y)
    
  
    
    if(missing(u)){
      u<-matrix(1,n,p)
    } else {      
      if(length(u)==p | length(u)==1 | identical(dim(u), c(n, p))){
        u <- matrix(u, n, p, byrow = is.vector(u))
        storage.mode(u) <- "double"
      }else {
        stop("Mispecified u, argument u must be either length 1 or p, or a n x p matrix,
          where p is the number of groups and n is the group size.")
      }
    }
    m<-nfactors+ncommon+p*ndistinct+p*nrandom
    
    
    P1<-diag(c(rep(1,nfactors),rep(0,ncommon+p*ndistinct),rep(1,p*nrandom)))
    if(!missing(random.cov)){
      if(!identical(dim(random.cov),c(nrandom,nrandom)))
        stop("Dimensions of the covariance matrix random.cov of random effects does not match with number of random effects.")
      P1[nfactors+ncommon+p*ndistinct+1:(p*nrandom),nfactors+ncommon+p*ndistinct+1:(p*nrandom)]<-
        as.matrix(.bdiag(replicate(p,random.cov,simplify=FALSE)))
      
    }
    
    #if(nrandom>0)
    #  P1[nfactors+ncommon+p*ndistinct+1:nrandom,nfactors+ncommon+p*ndistinct+1:nrandom]<-NA
    P1inf<-diag(c(rep(0,nfactors),rep(1,ncommon+p*ndistinct),rep(0,p*nrandom)))
    
    
    m1<-nfactors+ncommon+ndistinct+nrandom
    Z <- array(0,c(n,m1,p))
    Zind <- matrix(0,m1,p)
    tmp<-NULL
    for(i in 1:p){
      if(nfactors>0)
        tmp <- matrix(rep(factors[i,],n),ncol=nfactors,byrow=TRUE)
      Z[,,i]<-cbind(tmp,common_mm[,,i],distinct_mm[,,i],random_mm[,,i]) ##tänne nimet!!!
      Zind[,i]<-c(seq(from=1,to=nfactors+ncommon,length=nfactors+ncommon),
                  seq(from=nfactors+ncommon+1+(i-1)*ndistinct,to=nfactors+ncommon+i*ndistinct,length=ndistinct),
                  seq(from=nfactors+ncommon+p*ndistinct+1+(i-1)*nrandom,
                      to=nfactors+ncommon+p*ndistinct+i*nrandom,length=nrandom))
    }
    dimnames(Z)[[3]]<-names(splitdata)    
    dimnames(Z)[[2]]<-c(if(nfactors>0) paste("Factor",1:nfactors) else NULL,colnames(common_mm),colnames(distinct_mm),colnames(random_mm))
    
  
    model<-c(model,list(y=y,Z=aperm(Z,c(2,3,1)),Zind=Zind-1,a1=rep(0,m),P1=P1,P1inf=P1inf,
                        distribution=distribution,u=u,ncommon=ncommon,ndistinct=ndistinct,
                        nrandom=nrandom,nfactors=nfactors,tol=tol))
    
    class(model)<-"GLMMFmodel"
    model
  }