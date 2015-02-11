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

double approxF(const arma::mat& y, const arma::cube& Z, const arma::mat& u, const arma::vec& a1,
const arma::mat& P1, const arma::mat& P1inf, const int dist, const double tol, 
arma::mat&  ytilde, arma::mat& H, arma::mat& theta, const int maxiter, const int maxiter2, const double convtol, int& conv,
const arma::umat& zind, const int nfactors, const int trace) {
  
  int n = Z.n_slices;
  int p = Z.n_cols;
  arma::uvec ind(1);
  
  double ll_old;
  double ll=0.0;
  
  arma::mat coefs(a1.n_elem,Z.n_slices); 
  arma::mat coefs_old(a1.n_elem,Z.n_slices);
  arma::mat coefs_tmp(a1.n_elem,Z.n_slices);  
  switch (dist) {
    case 1:
    H = u;
    ytilde = y;
    break;
    case 2:
    H = exp(-theta)/u;
    ytilde = y%H + theta - 1.0;
    break;
    case 3:        
    H = (1.0+exp(theta))%(1.0+exp(theta))/(u%exp(theta));
    ytilde = theta + H*y - 1.0 - exp(theta);
    break;
    case 4:       
    H = 1.0/u;
    ytilde = theta+y/exp(theta)-1.0;
    break;
    case 5:
    H = (1.0/u+1.0/exp(theta));
    ytilde = theta+y/exp(theta)-1.0;
    break;
  }   
  
  ll_old = filterInApproxF(ytilde, Z, H, a1, P1, P1inf, tol,zind, nfactors, coefs);
  for(int t = 0; t<n; t++){
    ind =t;
    for(int i = 0; i<p; i++){
      theta(t,i) = arma::accu(Z.slice(t).col(i)%coefs(zind.col(i),ind));
    }
  }  
  
  if(dist > 1){
    ll_old += scaling(dist,y,u,theta,ytilde,H);
  } else return ll_old; //no need to iterate with gaussian model
  
  int iter2;
  int iter=0;
  while(iter < maxiter && arma::is_finite(ll)){  
    
    iter++;
    conv=iter;    
    switch (dist) {
      case 1:
      H = u;
      ytilde = y;
      break;
      case 2:
      H = exp(-theta)/u;
      ytilde = y%H + theta - 1.0;
      break;
      case 3:        
      H = (1.0+exp(theta))%(1.0+exp(theta))/(u%exp(theta));
      ytilde = theta + H*y - 1.0 - exp(theta);
      break;
      case 4:       
      H = 1.0/u;
      ytilde = theta+y/exp(theta)-1.0;
      break;
      case 5:
      H = (1.0/u+1.0/exp(theta));
      ytilde = theta+y/exp(theta)-1.0;
      break;
    }   
    
    ll = filterInApproxF(ytilde, Z, H, a1, P1, P1inf, tol,zind, nfactors, coefs);
    for(int t = 0; t<n; t++){
      ind =t;
      for(int i = 0; i<p; i++){
        theta(t,i) = arma::accu(Z.slice(t).col(i)%coefs(zind.col(i),ind));
      }
    } 
    if(dist > 1){
      ll += scaling(dist,y,u,theta,ytilde,H);
    }
    
    if( (((ll - ll_old)/(0.1 + std::abs(ll))) <= -convtol ) && iter>1 && maxiter2>0){
      iter2 = 0;      
      while(((ll - ll_old)/(0.1 + std::abs(ll))) < convtol && iter2<maxiter2 && arma::is_finite(ll)){
        iter2++;
        if(trace>0){
          Rcpp::Rcout<<"Step size halved due to decreasing likelihood."<<std::endl;
        }
        
        coefs = (coefs + coefs_old)/2.0;        
        for(int t = 0; t<n; t++){
          ind =t;
          for(int i = 0; i<p; i++){
            theta(t,i) = arma::accu(Z.slice(t).col(i)%coefs(zind.col(i),ind));
          }
        } 
        
        switch (dist) {
          case 1:
          H = u;
          ytilde = y;
          break;
          case 2:
          H = exp(-theta)/u;
          ytilde = y%H + theta - 1.0;
          break;
          case 3:        
          H = (1.0+exp(theta))%(1.0+exp(theta))/(u%exp(theta));
          ytilde = theta + H*y - 1.0 - exp(theta);
          break;
          case 4:       
          H = 1.0/u;
          ytilde = theta+y/exp(theta)-1.0;
          break;
          case 5:
          H = (1.0/u+1.0/exp(theta));
          ytilde = theta+y/exp(theta)-1.0;
          break;
        }
        ll = filterInApproxF(ytilde, Z, H, a1, P1, P1inf, tol,zind, nfactors,coefs_tmp);    
        for(int t = 0; t<n; t++){
          ind =t;
          for(int i = 0; i<p; i++){
            theta(t,i) = arma::accu(Z.slice(t).col(i)%coefs_tmp(zind.col(i),ind));
          }
        } 
        if(dist > 1){
          ll += scaling(dist,y,u,theta,ytilde,H);
        }        
      }
      if(iter2==maxiter2){
        if(trace>0){
          Rcpp::Rcout<<"Could not correct step size."<<std::endl;
        }       
      }      
      coefs = coefs_tmp;
      
    }         
    
    
    if(trace>1){
      Rcpp::Rcout<<"Iteration " << iter <<", Log-likelihood: " << ll <<std::endl;
    }
    
    
    
    
    
    if(!theta.is_finite()){
      conv = -3;
      break;      
    }
    if(std::abs(ll - ll_old)/(0.1 + std::abs(ll)) < convtol){
      break;
    }
    
    coefs_old = coefs;    
    ll_old = ll;
    
  }
  if(Rcpp::NumericMatrix::is_na(ll)){
    if(trace>0){
      Rcpp::Rcout << "Non-finite log-likelihood for Gaussian approximating model."<<std::endl;
    }
    conv = -2;
  }
  if(iter>=maxiter){
    if(trace>0){
      Rcpp::Rcout << "iteration limit reached in approximating algorithm."<<std::endl;
    }
    conv = 0;
  }
  
  return ll;
  
}

