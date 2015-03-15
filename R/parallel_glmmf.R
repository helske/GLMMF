#' Parallel estimation of GLMMF
#'
#' Function \code{parallel_glmmf} is a parallel version of glmmf function.
#'
#' @export
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
#' @param distribution Distribution of observations. Possible choices are "gaussian", "poisson", "binomial", "gamma and 
#' "negative binomial". Default is "gaussian".
#' @param init.random Initial values for random effect covariances.
#' @param init.dispersion Initial values for dispersion paremeters for Gaussian, negative binomial and Gamma distributions.
#' @param correlating.effects Logical. Default is TRUE.
#' @param maxiter Integer. Number of iterations for in iterative weighted least squares.

parallel_glmmf<-function(group,response,common.fixed,distinct.fixed,random,nfactors=0,data,
                         distribution=c("gaussian", "poisson", "binomial", "gamma", "negative binomial"),
                         u, correlating.effects=TRUE, common.dispersion=TRUE, init.random.cov, init.dispersion, init.factor,init.theta,
                         return.model=TRUE,maxiter=50,maxiter2=10, convtol=1e-8,estimate=TRUE,tol=.Machine$double.eps^0.5,trace=0,
                         seed=1,group_size=5,threads=1,samples=1,verbose=FALSE,...){  
  
  a<-proc.time()
  
  if(group_size<2)
    stop("Group size must be at least 2.")
  
  require(doSNOW) 
  set.seed(seed)
  
  
  mc<-match.call(expand.dots=TRUE)  
  mc$data<-NULL
  mc$init.dispersion<-NULL
  mc$init.factor<-NULL
  mc$init.theta<-NULL
  mc$seed<-NULL
  mc$group_size<-NULL
  mc$threads<-NULL
  mc$samples<-NULL
  mc$return.model<-FALSE
  mc$trace<-0
  mc$estimate<-TRUE
  mc[[1L]] <- quote(glmmf)
  
  
  fun <- function(i) {
    mc$data<-quote(data[data[,group]%in%groups[[i]],])
    mc$init.factor<-quote(init.factor[group_names%in%groups[[i]],])
    
    fit<-eval(mc)
    fit$fit$par
  }     
  
  group_names<-levels(data[,group])
  p<-length(group_names)
  n_groups<-p%/%group_size 
  
  
  
  cl <- makeCluster(threads)
  registerDoSNOW(cl)  
  initfactors_list<-foreach(iter=1:samples,.packages="foreach",.verbose=verbose,.export=c("mc","glmmf","group","data","init.factor")) %dopar%{    
    groups <- split(sample(group_names[-(1:(nfactors-1))]), cut(1:(p-nfactors+1), n_groups, labels = FALSE))  
    groups <- lapply(groups,append,group_names[1:(nfactors-1)],0)        
    result <- foreach(i=1:length(groups),.verbose=verbose) %do% fun(i)  
    initfactors<-matrix(0,p,nfactors)
    for(i in 1:length(groups)){
      Z<-matrix(0,nrow=length(groups[[i]]),ncol=nfactors)
      Z[lower.tri(Z,TRUE)]<-result[[i]]
      initfactors[group_names%in%(groups[[i]][-(1:(nfactors-1))]),]<-Z[-(1:(nfactors-1)),]  
      initfactors[1:(nfactors-1),]<-initfactors[1:(nfactors-1),]+Z[(1:(nfactors-1)),]
    }    
    initfactors[1:(nfactors-1),]<-initfactors[1:(nfactors-1),]/length(groups)
    initfactors[upper.tri(initfactors)]<-0
    initfactors
  } 
  
  print("Submodelling completed.")
  mc$data<-quote(data)
  mc$estimate=FALSE
  
  fun <- function(i) {
    mc$init.factor<-initfactors_list[[i]]
    fit<-eval(mc)
    #initfactors<-matrix(0,p,nfactors)
    #initfactors[lower.tri(initfactors,TRUE)]<-fit$fit$par
    #list(aic=-2.0*fit$logLik+2.0*(nfactors*p-1),initfactors=initfactors)    
    fit$logLik
  }   
  
  res<-foreach(iter=1:samples,.export=c("glmmf","mc","data"),.combine="c",.verbose=verbose) %dopar% fun(iter)
  stopCluster(cl)
  
  print("Parallel estimation of samples completed.")
  
  #aic <- sapply(res, "[[", 1)
  #w<- exp(-0.5*(aic - min(aic)))
  #w<-w/sum(w)
  #initfactors<-Reduce("+", Map("*",lapply(res, "[[", 2), w))
  initfactors <- initfactors_list[[which.max(res)]]
  f1<-glmmf(group=group,response=response,common.fixed=common.fixed,distinct.fixed=distinct.fixed,random=random,
           nfactors=nfactors,data=data,
           distribution=distribution,
           u=u, correlating.effects=correlating.effects, common.dispersion=common.dispersion, init.random.cov=init.random.cov, 
           init.dispersion=init.dispersion, init.factor=initfactors,init.theta=init.theta,
           return.model=return.model,maxiter=maxiter,maxiter2=maxiter2, convtol=convtol,estimate=FALSE,tol=tol,trace=1,...) 
  
  print(f1$logLik)
  print(proc.time()-a)
  
  f<-glmmf(group=group,response=response,common.fixed=common.fixed,distinct.fixed=distinct.fixed,random=random,
        nfactors=nfactors,data=data,
        distribution=distribution,
        u=u, correlating.effects=correlating.effects, common.dispersion=common.dispersion, init.random.cov=init.random.cov, 
        init.dispersion=init.dispersion, init.factor=initfactors,init.theta=f$linear.predictor,
        return.model=return.model,maxiter=maxiter,maxiter2=maxiter2, convtol=convtol,estimate=estimate,tol=tol,trace=trace,...) #initmat=-1499, optim=-1346.19
  f
}
