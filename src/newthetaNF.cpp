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

double newthetaNF(const arma::mat& y, const arma::cube& Z, const arma::mat& H,
const arma::vec& a1, const arma::mat& P1, const arma::mat& P1inf, const double tol,
const arma::umat& zind, arma::mat& theta) {
  
  

  double tiny = std::numeric_limits<double>::min();
  
  arma::mat pt(P1.begin(),P1.n_rows,P1.n_rows);  
  at = a1;
  arma::vec kt(pt.n_rows);
  
  
  unsigned int m1 = Z.n_rows;
  
  double vt, ft;
  unsigned int j=0;
  unsigned int d=0;
  int rankp = arma::accu(P1inf);
  
  double lik = 0.0;
  const double c = 0.5*std::log(2.0*M_PI);
  
  kt = pt.cols(zind.col(0))*Z.slice(0).col(0);  
  double yhat = arma::accu(Z.slice(0).col(0)%at.rows(zind.col(0)));
  double zk = arma::accu(Z.slice(0).col(0)%kt.rows(zind.col(0)));
  
  
  
  if(rankp>0){
    
    arma::vec kinf(pt.n_rows);
    arma::mat pinf(P1inf.begin(),pt.n_rows,pt.n_rows);
    double finf;
    
    kinf = pinf.cols(zind.col(0))*Z.slice(0).col(0);   
    finf = arma::accu(Z.slice(0).col(0)%kinf.rows(zind.col(0)));  
    while(d < Z.n_slices && rankp>0){   
      
      for(j=0; j<Z.n_cols; j++){
        if(arma::is_finite(y(d,j))){
          ft = zk + H(d,j);
          if(finf>tol){
            finf = 1.0/finf;
            lik += 0.5*log(finf);
            vt = y(d,j)-yhat;
            at = at + kinf*vt*finf;     
            //pinf = pinf - kinf*kinf.t()*finf;
            //pt = pt +(kinf*kinf.t()*ft*finf-kt*kinf.t()-kinf*kt.t())*finf; 
            for(unsigned int k = 0; k<pinf.n_rows; k++){
              for(unsigned int l = k; l<pinf.n_rows; l++){                
                pinf(l,k) = pinf(l,k) - kinf(l)*kinf(k)*finf; 
                pt(l,k) = pt(l,k) +kinf(l)*kinf(k)*pow(finf,2)*ft-
                kt(l)*kinf(k)*finf-kinf(l)*kt(k)*finf;
              }
            }
            //pt = symmatl(pt);
            //pinf = symmatl(pinf);  
            
            rankp--;
          } else {
            if (ft > tiny){
              ft=1.0/ft;
              vt = y(d,j)-yhat;
              lik -= c + 0.5*(-log(ft) + pow(vt,2)*ft);
              at = at + kt*vt*ft;
              //pt = pt - kt*kt.t()*ft; 
              for(unsigned int k = 0; k<pt.n_rows; k++){
                for(unsigned int l = k; l<pt.n_rows; l++){
                  pt(l,k) = pt(l,k) - kt(l)*kt(k)*ft;          
                }
              }
              //pt = symmatl(pt); 
            }
          }
        }
        if(j==(Z.n_cols-1)){

          if(d<(Z.n_slices-1)){
            //kt =     pt.cols(zind.col(0))*Z.slice(d+1).col(0);   
            //kinf = pinf.cols(zind.col(0))*Z.slice(d+1).col(0);
            kt.zeros();
            for(unsigned int k = 0; k<pt.n_rows; k++){
              for(unsigned int l = 0; l<m1; l++){
                kt(k) = kt(k) + pt(std::max(k,zind(l,0)),std::min(k,zind(l,0)))*Z(l,0,d+1);       
              }
            }  
            kinf.zeros();
            for(unsigned int k = 0; k<pt.n_rows; k++){
              for(unsigned int l = 0; l<m1; l++){
                kinf(k) = kinf(k) + pinf(std::max(k,zind(l,0)),std::min(k,zind(l,0)))*Z(l,0,d+1);       
              }
            }  
            yhat = arma::accu(Z.slice(d+1).col(0)%at.rows(zind.col(0)));
            zk = arma::accu(Z.slice(d+1).col(0)%kt.rows(zind.col(0)));          
            finf = arma::accu(Z.slice(d+1).col(0)%kinf.rows(zind.col(0)));  
          }
        } else {
          //next observation
          //kt =  pt.cols(zind.col(j+1))*Z.slice(d).col(j+1);
          //kinf = pinf.cols(zind.col(j+1))*Z.slice(d).col(j+1);
          kt.zeros();
          for(unsigned int k = 0; k<pt.n_rows; k++){
            for(unsigned int l = 0; l<m1; l++){
              kt(k) = kt(k) + pt(std::max(k,zind(l,j+1)), std::min(k,zind(l,j+1))) * Z(l,j+1,d);       
            }
          }
          kinf.zeros();
          for(unsigned int k = 0; k<pt.n_rows; k++){
            for(unsigned int l = 0; l<m1; l++){
              kinf(k) = kinf(k) + pinf(std::max(k,zind(l,j+1)), std::min(k,zind(l,j+1))) * Z(l,j+1,d);       
            }
          }
          yhat = arma::accu(Z.slice(d).col(j+1)%at.rows(zind.col(j+1)));
          zk = arma::accu(Z.slice(d).col(j+1)%kt.rows(zind.col(j+1))); 
          finf = arma::accu(Z.slice(d).col(j+1)%kinf.rows(zind.col(j+1)));
          
        }
        if(rankp==0){
          d--; // negate the last increment below 
          break;
        }
        
      }
      d++;
    }
    
    for(unsigned int i = j+1; i<Z.n_cols; i++){
      
      if(arma::is_finite(y(d,i))){
        ft = zk + H(d,i);
        
        if (ft > tiny){
          ft = 1.0/ft;
          vt = y(d,i)-yhat;
          lik -= c + 0.5*(-log(ft) + pow(vt,2)*ft);
          at = at + kt*vt*ft;
          for(unsigned int k = 0; k<pt.n_rows; k++){
            for(unsigned int l = k; l<pt.n_rows; l++){
              pt(l,k) = pt(l,k) - kt(l)*kt(k)*ft;          
            }
          }
          //pt = symmatl(pt);         
        }
      }
      
      if(i==(Z.n_cols-1)){
        //next time step
   
        if(d<(Z.n_slices-1)){
          
          //kt = pt.cols(zind.col(0))*Z.slice(d+1).col(0);     
          kt.zeros();
          for(unsigned int k = 0; k<pt.n_rows; k++){
            for(unsigned int l = 0; l<m1; l++){
              kt(k) = kt(k) + pt(std::max(k,zind(l,0)),std::min(k,zind(l,0)))*Z(l,0,d+1);       
            }
          }   
          yhat = arma::accu(Z.slice(d+1).col(0)%at.rows(zind.col(0)));
          zk = arma::accu(Z.slice(d+1).col(0)%kt.rows(zind.col(0)));         
        }
      } else {
        //next observation
        // kt = pt.cols(zind.col(i+1))*Z.slice(d).col(i+1);      
        kt.zeros();
        for(unsigned int k = 0; k<pt.n_rows; k++){
          for(unsigned int l = 0; l<m1; l++){
            kt(k) = kt(k) + pt(std::max(k,zind(l,i+1)), std::min(k,zind(l,i+1))) * Z(l,i+1,d);       
          }
        }
        yhat = arma::accu(Z.slice(d).col(i+1)%at.rows(zind.col(i+1)));
        zk = arma::accu(Z.slice(d).col(i+1)%kt.rows(zind.col(i+1))); 
        
      }
      
      
    }
    d++;
  }    
  
  
  
  for(unsigned int t = d; t<Z.n_slices; t++){   
    for(unsigned int i = 0; i<Z.n_cols; i++){
      if(arma::is_finite(y(t,i))){
        ft = zk + H(t,i);
        
        if (ft > tiny){
          ft = 1.0/ft;
          vt = y(t,i)-yhat;
          lik -= c + 0.5*(-log(ft) + pow(vt,2)*ft);
          at = at + kt*vt*ft;
          for(unsigned int k = 0; k<pt.n_rows; k++){
            for(unsigned int l = k; l<pt.n_rows; l++){
              pt(l,k) = pt(l,k) - kt(l)*kt(k)*ft;          
            }
          }         
        }
      }
      if(i==(Z.n_cols-1)){
        //next time step
        if(t<(Z.n_slices-1)){  
          kt.zeros();
          for(unsigned int k = 0; k<pt.n_rows; k++){
            for(unsigned int l = 0; l<m1; l++){
              kt(k) = kt(k) + pt(std::max(k,zind(l,0)),std::min(k,zind(l,0)))*Z(l,0,t+1);       
            }
          }        
          yhat = arma::accu(Z.slice(t+1).col(0)%at.rows(zind.col(0)));
          zk = arma::accu(Z.slice(t+1).col(0)%kt.rows(zind.col(0)));         
        }
      } else {
        //next observation      
        kt.zeros();
        for(unsigned int k = 0; k<pt.n_rows; k++){
          for(unsigned int l = 0; l<m1; l++){
            kt(k) = kt(k) + pt(std::max(k,zind(l,i+1)), std::min(k,zind(l,i+1))) * Z(l,i+1,t);       
          }
        }  
        yhat = arma::accu(Z.slice(t).col(i+1)%at.rows(zind.col(i+1)));
        zk = arma::accu(Z.slice(t).col(i+1)%kt.rows(zind.col(i+1))); 
        
      }        
      
    } 
  }
  
  for(int t = 0; t<n; t++){
    ind =t;
    for(int i = 0; i<p; i++){
      theta(t,i) = arma::accu(Z.slice(t).col(i)%at(zind.col(i),ind));
    }
  }  
  
  return lik;
}