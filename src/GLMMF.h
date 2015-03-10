#ifndef GLMMFINTERNALS
#define GLMMFINTERNALS

//#define ARMA_NO_DEBUG
#include "RcppArmadillo.h"

void smootherF(const arma::mat& y, const arma::cube& Z, const arma::mat& H,
const arma::colvec& a1, const arma::mat& P1, const arma::mat& P1inf, const double tol,
const arma::umat& zind, const int nfactors, arma::mat& coefs,arma::cube& coefVars);

void smootherNF(const arma::mat& y, const arma::cube& Z, const arma::mat& H,
const arma::colvec& a1, const arma::mat& P1, const arma::mat& P1inf, const double tol,
const arma::umat& zind, arma::mat& coefs,arma::cube& coefVars);


double newthetaF(const arma::mat& y, const arma::cube& Z, const arma::mat& H,
const arma::vec& a1, const arma::mat& P1, const arma::mat& P1inf, const double tol, 
const arma::umat& zind, const int nfactors, arma::mat& theta);

double newthetaNF(const arma::mat& y, const arma::cube& Z, const arma::mat& H,
const arma::vec& a1, const arma::mat& P1, const arma::mat& P1inf, const double tol,
const arma::umat& zind, arma::mat& theta);

double approxF(const arma::mat& y, const arma::cube& Z, const arma::mat& u, const arma::vec& a1, 
const arma::mat& P1, const arma::mat& P1inf, const int dist, const double tol, 
arma::mat&  ytilde, arma::mat& H, arma::mat& theta, const int maxiter, const int maxiter2, const double convtol, int& conv,
const arma::umat& zind, const int nfactors, const int trace);

double approxNF(const arma::mat& y, const arma::cube& Z, const arma::mat& u, const arma::vec& a1, 
const arma::mat& P1, const arma::mat& P1inf, const int dist, const double tol, 
arma::mat&  ytilde, arma::mat& H, arma::mat& theta, const int maxiter, const int maxiter2, const double convtol, int& conv,
const arma::umat& zind, const int trace);

double scaling(const int dist, const arma::mat& y, const arma::mat& u,const arma::mat& theta,
               const arma::mat& ytilde, const arma::mat& H);

void ytildeH(const int dist, const arma::mat& y, const arma::mat& u, const arma::mat& theta, arma::mat& ytilde, arma::mat& H);

double pytheta(const int dist, const arma::mat& y, const arma::mat& u,const arma::mat& theta);
               
double ptheta(const arma::mat& y, const arma::cube& Z, const arma::mat& H,
const arma::colvec& a1, const arma::mat& P1, const arma::mat& P1inf, const double tol,
const arma::umat& zind, const int nfactors);

#endif