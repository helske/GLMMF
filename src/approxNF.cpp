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

double approxNF(const arma::mat& y, const arma::cube& Z, const arma::mat& u, const arma::vec& a1,
const arma::mat& P1, const arma::mat& P1inf, const int dist, const double tol, 
arma::mat&  ytilde, arma::mat& H, arma::mat& theta, const int maxiter, const int maxiter2, const double convtol, int& conv,
const arma::umat& zind, const int trace) {
  
  
  
  double ll_old;
  double ll=0.0;
  
  arma::vec at(a1.n_elem);
  arma::vec at_old(a1.n_elem);
  arma::vec at_tmp(a1.n_elem);
  
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
  
  
  
  ll_old = filterInApproxNF(ytilde, Z, H, a1, P1, P1inf, tol, zind, at);
  
  for(unsigned int t = 0; t<Z.n_slices; t++){
    for(unsigned int i = 0; i<Z.n_cols; i++){
      theta(t,i) = arma::accu(Z.slice(t).col(i)%at.rows(zind.col(i)));     
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
    
    
    
    ll = filterInApproxNF(ytilde, Z, H, a1, P1, P1inf, tol, zind, at);  
    for(unsigned int t = 0; t<Z.n_slices; t++){
      for(unsigned int i = 0; i<Z.n_cols; i++){
        theta(t,i) = arma::accu(Z.slice(t).col(i)%at.rows(zind.col(i)));     
      }
    }
    if(dist > 1){
      ll += scaling(dist,y,u,theta,ytilde,H);
    }   
    
    if( (((ll - ll_old)/(0.1 + std::abs(ll))) <= -convtol ) && iter>1 && maxiter2>0 && arma::is_finite(ll)) {
      iter2 = 0;      
      while(((ll - ll_old)/(0.1 + std::abs(ll))) < convtol && iter2<maxiter2){
        iter2++;
        if(trace>0){
          Rcpp::Rcout<<"Step size halved due to decreasing likelihood."<<std::endl;
        }       
        at = (at + at_old)/2.0;
        for(unsigned int t = 0; t<Z.n_slices; t++){
          for(unsigned int i = 0; i<Z.n_cols; i++){
            theta(t,i) = arma::accu(Z.slice(t).col(i)%at.rows(zind.col(i)));     
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
        ll = filterInApproxNF(ytilde, Z, H, a1, P1, P1inf, tol, zind, at_tmp);    
        for(unsigned int t = 0; t<Z.n_slices; t++){
          for(unsigned int i = 0; i<Z.n_cols; i++){
            theta(t,i) = arma::accu(Z.slice(t).col(i)%at_tmp.rows(zind.col(i)));     
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
    at_old = at;    
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
