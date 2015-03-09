// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "GLMMF.h"
// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// simple example of creating two matrices and
// returning the result of an operatioon on them
//
// via the exports attribute we tell Rcpp to make this function
// available from R
// [[Rcpp::export]]

Rcpp::List logLik(const Rcpp::NumericMatrix y_, const Rcpp::NumericVector Z_, const Rcpp::NumericMatrix u_,
const Rcpp::NumericVector a1_, const Rcpp::NumericMatrix P1_, const Rcpp::NumericMatrix P1inf_, 
const int dist, double tol, int maxiter, int maxiter2, double convtol, Rcpp::NumericMatrix theta_
, const Rcpp::IntegerMatrix Zind_, const int nfactors, const int trace,const int compgrad){
  
  const unsigned int n = y_.nrow();
  const unsigned int p = y_.ncol();
  const unsigned int m = P1_.nrow();
  const unsigned int m1 = Zind_.nrow();
  
  const arma::mat y(y_.begin(), n,p, false);
  const arma::cube Z(Z_.begin(), m1, p, n, false);
  const arma::mat u(u_.begin(), n,p, false);
  const arma::vec a1(a1_.begin(),m, false);
  const arma::mat P1(P1_.begin(),m,m, false);
  const arma::mat P1inf(P1inf_.begin(),m,m, false);
  
  const arma::imat izind(Zind_.begin(),m1,p,false);
  const arma::umat zind = arma::conv_to<arma::umat>::from(izind);
  
  arma::mat theta(theta_.begin(), n,p, true);
  arma::mat ytilde(n,p);
  arma::mat H(n,p);
  
  double loglik; 
  int conv = 0; 
  
  // approximating Gaussian model 
  if(nfactors>0){
    loglik = approxF(y, Z, u, a1, P1, P1inf, dist, tol, ytilde, H, theta, maxiter, maxiter2, convtol, 
    conv, zind, nfactors,trace);   
  } else {  
    loglik = approxNF(y, Z, u, a1, P1, P1inf, dist, tol, ytilde, H, theta, maxiter, maxiter2, convtol, 
    conv, zind, trace);
  }
  if(conv<0){
    return Rcpp::List::create(Rcpp::Named("logLik") = -std::numeric_limits<double>::max(), Rcpp::Named("convergence") = conv);
  }  
  
  logLik =+ scaling(dist, y, u, theta, ytilde, H);
  
  arma::mat grad(p,nfactors);
  grad.zeros();
  if(compgrad==1 && nfactors>0){
    arma::mat coefs(m,n);
    arma::cube coefVars(m,m,n); 
    arma::rowvec tmpm(m);
    arma::rowvec tmpm2(m);
    if(nfactors>0){
      smootherF(ytilde, Z, H, a1, P1, P1inf, tol, zind, nfactors,coefs,coefVars);
    } else {
      smootherNF(ytilde, Z, H, a1, P1, P1inf, tol, zind, coefs, coefVars);
    }
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
    
  }
  
  
  
  return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("gradient") = grad,
  Rcpp::Named("logLik") = loglik,Rcpp::Named("convergence") = conv);
}
