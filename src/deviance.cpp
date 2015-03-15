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

double deviance(const arma::mat& y, const arma::mat& u, 
const Rcpp::IntegerVector& dist, const arma::mat& theta) {
  
  
  // This function is not optimized at all!!!
  
  const unsigned int n = y.n_rows;
  const unsigned int p = y.n_cols; 
  arma::mat res(n,p);  
  arma::mat mu(n,p);
  arma::vec tmp(n);
  arma::vec tmp2(n);
  
            
  res = y;
  
  for(int i = 0; i<dist.size(); i++){      
    switch (dist(i)) {
      case 1:
      res.col(i) = arma::pow(res.col(i)-theta.col(i),2);      
      break;
      case 2:
      mu.col(i) =  exp(theta.col(i))%u.col(i);    
     
 
      tmp = res.col(i)/mu.col(i);
      tmp.elem(find(res.col(i) == 0.0) ).ones();
      res.col(i) =  2.0*(res.col(i)%log(tmp) - res.col(i) + mu.col(i));      
      
//      // from glm:  good <- (weights > 0) & (mu.eta.val != 0)
//Rcpp::Rcout<<res.col(i)<<std::endl;
//      tmp = res.col(i);      
//      tmp.elem(find(mu.col(i)<1e-12)).zeros();
//      res.col(i) = tmp;
     
      break;
      case 3:        
      mu.col(i) = exp(theta.col(i))/(1.0+exp(theta.col(i)));
      
      res.col(i) = res.col(i)/u.col(i);
      tmp = res.col(i)/mu.col(i);
      tmp.elem(find(res.col(i) == 0.0) ).ones();                
      
      tmp2 = (1.0 - res.col(i))/(1 - mu.col(i));
      tmp2.elem(find(res.col(i) == 1.0) ).ones(); 
      tmp2.elem(find(mu.col(i) == 1.0) ).ones();                
      res.col(i) = 2.0 * u.col(i) % (res.col(i) % log(tmp) + (1.0 - 
      res.col(i)) % log(tmp2));
      break;
      case 4:
      mu.col(i) = exp(theta.col(i));
      
      tmp = res.col(i)/mu.col(i);
      tmp.elem(find(res.col(i) == 0.0) ).ones();
      
      res.col(i) = -2.0 * (log(tmp) - (res.col(i)-mu.col(i))/mu.col(i));      
      
      break;
      case 5:
      mu.col(i) = exp(theta.col(i));
      
      tmp = res.col(i)/mu.col(i);
      for(unsigned int j=1;j<n;j++){
        if(res(j,i)<1.0){
          tmp(j,i) = 1.0/mu(j,i);
        }
      }
      
      res.col(i) = 2.0 * (res.col(i)%log(tmp) - 
      (res.col(i)+u.col(i))%log((res.col(i)+u.col(i))/(mu.col(i)+u.col(i))));
      break;
    }
  }  
  return arma::accu( res.elem( find_finite(res) ) );
  
  
}