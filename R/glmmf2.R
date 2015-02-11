#' it Generalized Linear Mixed-Effects Models with Latent Factors Using State
#' Space Framework with Approximate Speedup
#' 
#' Function \code{glmmf2} estimates GLMM with latent factors using methods based
#' on state space modelling with speed up method described in [1].
#' 
#' @export
#' @import MASS
#' @importFrom Rcpp evalCpp
#' @useDynLib GLMMF
#' @param group Name of the grouping variable in data. Only one grouping
#'   variable is allowed and the group sizes must be equal. In case of unequal
#'   group sizes, patch missing rows with NA's.
#' @param response Name of the response variable in data.
#' @param common.fixed formula for common fixed effects. LHS of the formula is
#'   ignored if present.
#' @param distinct.fixed Formula for distinct fixed effects i.e. each group has
#'   separate regression coefficient. This formula cannot contain variables
#'   which are already present in \code{common.fixed} or \code{random} formulas,
#'   as in that case the model would not be identifiable. LHS of the formula is
#'   ignored if present.
#' @param random Formula for random effects. LHS of the formula is ignored if
#'   present.
#' @param data Data frame containing the variables in the model. Must contain
#'   all variables used formulas and variables defined in \code{group} and
#'   \code{response}
#' @param distribution Distribution of observations. Possible choices are
#'   "gaussian", "poisson", "binomial", "gamma and "negative binomial". Default
#'   is "gaussian".
#' @param init.random Initial values for random effect covariances.
#' @param init.dispersion Initial values for dispersion paremeters for Gaussian,
#'   negative binomial and Gamma distributions.
#' @param correlating.effects Logical. Default is TRUE.
#' @param nsim Integer. Number of independent samples used in importance
#'   sampling. Default is 0, which corresponds to Laplace approximation. Not yet
#'   implemented.
#' @param maxiter Integer. Number of iterations for in iterative weighted least
#'   squares.

glmmf2<-
  function(group,response,common.fixed,distinct.fixed,random,nfactors=0,data,
           distribution=c("gaussian", "poisson", "binomial", "gamma", "negative binomial"),
           u, correlating.effects=TRUE, common.dispersion=TRUE, init.random.cov, init.dispersion, init.factor,init.theta=NULL,
           nsim=0, return.model=TRUE,maxiter=50,maxiter2=10, convtol=1e-8,estimate=TRUE,
           tol=.Machine$double.eps^0.5,trace=0,gradient=TRUE,maxapp=10,epsilon=1e-5,...){
    
  
    distribution<-match.arg(distribution)
    estimate.dispersion <- distribution%in%c("gaussian","negative binomial","gamma")
    
    if(missing(init.dispersion) && estimate.dispersion){      
     
      init.dispersion<-
        switch(distribution,
               "gaussian" = {
                 if(!missing(common.fixed)){
                   tmpf<-formula(paste(response,paste(common.fixed,collapse="")))
                 } else tmpf<-formula(paste(response,"~1"))       
                 summary(glm(tmpf,data=data,family=gaussian))$dispersion
               },"negative binomial" = {
                 if(!missing(common.fixed)){
                   tmpf<-formula(paste(response,paste(common.fixed,collapse="")))
                 } else tmpf<-formula(paste(response,"~1"))       
                 glm.nb(tmpf,data=data)$theta
               },"gamma"= {
                 if(!missing(common.fixed)){
                   tmpf<-formula(paste(response,paste(common.fixed,collapse="")))
                 } else tmpf<-formula(paste(response,"~1"))       
                 summary(glm(tmpf,data=data,family=Gamma(link="log")))$dispersion
               })
     
      
    }
    
    model<-buildGLMMF(group=group,response=response,common.fixed=common.fixed,
                      distinct.fixed=distinct.fixed,random=random,nfactors=nfactors,data=data,
                      u=u,distribution=distribution,tol=tol,random.cov=init.random.cov,factors=init.factor)
    
    
    
    res<-fitGLMMF2(model=model,estimate.dispersion=estimate.dispersion,common.dispersion=common.dispersion,
                  correlating.effects=correlating.effects,maxiter=maxiter,maxiter2=maxiter2, convtol=convtol,
                  init.random.cov=init.random.cov,init.dispersion=init.dispersion,
                  init.factor=init.factor,init.theta=init.theta,estimate=estimate,trace,gradient=gradient,maxapp=maxapp,epsilon=epsilon,...)    
    model<-res$model
    
    if(is.null(init.theta) || init.theta=="iterative"){
    theta <- init_theta(model$y, model$u, model$distribution) 
    } else theta <-init.theta
    dist<-pmatch(x = model$distribution, 
                     table = c("gaussian", "poisson", "binomial", 
                               "gamma", "negative binomial"))
    out<-filterNoSim(model$y, model$Z, model$u, model$a1, model$P1, model$P1inf, dist, 
                     model$tol, maxiter, maxiter2, convtol, theta, model$Zind, model$nfactors,trace)
   
    
    if(out$conv<0){
      stop(paste("Approximating algorithm did not converge. Error code ",out$conv))      
    }  results<-c(out,res)
    class(results)<-"glmmf.results"
    results
    
    
    
  }