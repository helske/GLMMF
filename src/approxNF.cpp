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
arma::mat&  ytilde, arma::mat& H, arma::mat& theta, const int maxiter, const int maxiter2, 
const double convtol, int& conv,
const arma::umat& zind, const int trace) {
  
  int n = Z.n_slices;
  int p = Z.n_cols;
  arma::uvec ind(1);
  double lik;
  
  
  double dev_old = 0.0;
  double dev;
  
    
  dev_old = pytheta(dist,y,u,theta,ytilde,H)+ ptheta(theta, Z, a1, P1, P1inf, tol, zind, 0);
 
  if(dist > 1){
    return dev_old; //no need to iterate with gaussian model
  }
  
  arma::mat theta_old(n,p);    
  theta_old = theta;
  
  
  
  
  int iter2;
  int iter  =0;
  while(iter < maxiter){  
    
    iter++;
    conv = iter;   
    
    ytildeH(dist, y, u, theta, ytilde, H);
    lik = newthetaNF(ytilde, Z, H, a1, P1, P1inf, tol, zind, theta);
    dev = pytheta(dist,y,u,theta,ytilde,H) + ptheta(theta, Z, a1, P1, P1inf, tol, zind, 0);
    
    if( (((dev - dev_old)/(0.1 + std::abs(dev))) < convtol ) && iter>1 && maxiter2>0){
      iter2 = 0;      
      while(((dev - dev_old)/(0.1 + std::abs(dev))) < 0.0 && iter2<maxiter2){
        iter2++;
        if(trace>0){
          Rcpp::Rcout<<"Step size halved due to decreasing likelihood."<<std::endl;
        }
        
        theta = 0.5*(theta + theta_old);        
        
        ytildeH(dist, y, u, theta, ytilde, H);
        lik = newthetaF(ytilde, Z, H, a1, P1, P1inf, tol, zind, nfactors, theta);
        dev = pytheta(dist,y,u,theta,ytilde,H) + ptheta(theta, Z, a1, P1, P1inf, tol, zind, 0);
        
      }
      if(iter2==maxiter2){
        if(trace>0){
          Rcpp::Rcout<<"Could not correct step size."<<std::endl;
        }       
      }      
      
    }         
    
    if(trace>1){
      Rcpp::Rcout<<"Iteration " << iter <<", p(theta|y): " << dev <<std::endl;
    }
    
    if(!theta.is_finite()){
      conv = -3;
      break;      
    }
    if(std::abs(dev - dev_old)/(0.1 + std::abs(dev)) < convtol){
      break;
    } else {
      theta_old = theta;
      dev_old = dev;
    }
    
    
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
  
  return lik;
}