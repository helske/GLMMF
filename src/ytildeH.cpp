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
void ytildeH(const int dist, const arma::mat& y, const arma::mat& u, const arma::mat& theta, arma::mat& ytilde, arma::mat& H) {
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
}
