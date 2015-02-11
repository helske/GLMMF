p<-12
n<-28

gradloglik<-function(pars,model){
  p<-attr(model,"p")
  n<-attr(model,"n")
  out<-KFS(model,filtering="state",smoothing="state")
  gradient<-matrix(0,p,2)
  
  for(j in 1:p){
    for(k in 1:2){
      if(j>=k){
        for(t in 1:n){
          gradient[j,k] <- gradient[j,k] + 
            (model$y[t,j]*out$alpha[t,12+k]-(out$V[12+k,,t]+out$alpha[t,12+k]%*%t(out$alpha[t,]))%*%model$Z[j,,t])/model$H[j,j,t]
        }
      }
    }
  }
  ll<- -out$logLik
  attr(ll,"gradient")<- -gradient[lower.tri(gradient,TRUE)]
  ll
}


model<-approxSSM(model)


library(numDeriv)
likfn<-function(pars){
  Z<-matrix(0,nrow=p,ncol=2)
  Z[lower.tri(Z,TRUE)]<-pars
  model$Z[,13:14,]<-Z
  logLik(model)
}
grad(likfn,model$Z[,13:14,1][lower.tri(matrix(0,p,2),TRUE)])

for(int j=0; j<p; j++){
  for(int k=0; k<nfactors; k++){
    if(j>=k){
      for(int t=0; t<n; t++){
        tmpm = coefVars.slice(t).row(k);
        tmpm2 = coefs.col(t).t();
        
        grad(j,k) += (ytilde(t,j)*coefs(k,t)-
                        arma::as_scalar((tmpm.cols(zind.col(j)) + coefs(k,t)*tmpm2.cols(zind.col(j)))*Z.slice(t).col(j)))/H(t,j);              
      }
    }
  }
}