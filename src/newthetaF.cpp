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

double newthetaF(const arma::mat& y, const arma::cube& Z, const arma::mat& H,
const arma::colvec& a1, const arma::mat& P1, const arma::mat& P1inf, const double tol,
const arma::umat& zind, const int nfactors, arma::mat& theta) {
  
  
  double tiny = pow(std::numeric_limits<double>::epsilon(),0.75);
  
  int n = Z.n_slices;
  int p = Z.n_cols;
  unsigned int m = P1.n_rows;
  unsigned int m1 = Z.n_rows;  
  int j=0;
  int d=-1;
  arma::vec at(a1.begin(),m);           // m
  //at = a1;
  arma::mat coefs(a1.n_elem,Z.n_slices);  
  arma::mat vt(p,n);            // p x n
  
  arma::mat pt(P1.begin(),m,m); // m x m
  arma::cube kt(m,p,n); // m x p x n
  arma::mat ft(p,n);            // p x n
  
  
  arma::mat pinf(P1inf.begin(),m,m); // m x m
  arma::cube kinf(m,p,n); // m x p x n
  arma::mat finf(p,n);            // p x n
  
  arma::cube PtAll(m,m,n);
  arma::cube PinfAll(m,m,n);
  
  coefs.col(0)=a1;
  PinfAll.slice(0) = P1inf; 
  PtAll.slice(0) = P1;   
  
  int rankp = arma::accu(P1inf);
  
  double lik = 0.0;
  const double c = 0.5*std::log(2.0*M_PI);
  
  double finv;
  arma::vec tmpm(m);
  
  
  kt.zeros();
  kinf.zeros();
  kt.slice(0).col(0) = pt.cols(zind.col(0))*Z.slice(0).col(0);    
  
  double yhat = arma::accu(Z.slice(0).col(0)%at.rows(zind.col(0)));
  tmpm = kt.slice(0).col(0);
  double zk = arma::accu(Z.slice(0).col(0)%tmpm.rows(zind.col(0)));  
  
  if(rankp>0){
    
    kinf.slice(0).col(0) = pinf.cols(zind.col(0))*Z.slice(0).col(0); 
    tmpm = kinf.slice(0).col(0);
    finf(0,0) = arma::accu(Z.slice(0).col(0)%tmpm.rows(zind.col(0)));  
    while(d < n && rankp>0){   
      d++;
      for(j=0; j<p; j++){
        if(arma::is_finite(y(d,j))){
          ft(j,d) = zk + H(d,j);
          if(finf(j,d)>tol){
            finv = 1.0/finf(j,d);
            lik += 0.5*log(finv);
            vt(j,d) = y(d,j)-yhat;
            at = at + kinf.slice(d).col(j)*vt(j,d)*finv;             
            //pinf = pinf - kinf*kinf.t()*finf;
            //pt = pt +(kinf*kinf.t()*ft*finf-kt*kinf.t()-kinf*kt.t())*finf; 
            for(unsigned int k = 0; k<m; k++){
              for(unsigned int l = k; l<m; l++){                
                pinf(l,k) = pinf(l,k) - kinf(l,j,d)*kinf(k,j,d)*finv; 
                pt(l,k) = pt(l,k) +kinf(l,j,d)*kinf(k,j,d)*pow(finv,2)*ft(j,d)-
                kt(l,j,d)*kinf(k,j,d)*finv-kinf(l,j,d)*kt(k,j,d)*finv;
              }
            }
            //pt = symmatl(pt);
            //pinf = symmatl(pinf);  
            
            rankp--;
          } else {
            if (ft(j,d) > tiny){
              finv=1.0/ft(j,d);
              vt(j,d) = y(d,j)-yhat;
              lik -= c + 0.5*(-log(finv) + pow(vt(j,d),2)*finv);
              at = at + kt.slice(d).col(j)*vt(j,d)*finv;
              //pt = pt - kt*kt.t()*ft; 
              for(unsigned int k = 0; k<m; k++){
                for(unsigned int l = k; l<m; l++){
                  pt(l,k) = pt(l,k) - kt(l,j,d)*kt(k,j,d)*finv;         
                }
              }
              //pt = symmatl(pt); 
            }
          }
        }
        //        
        if(j==(p-1)){
          //next time step
          if(d<(n-1)){
            
            pt.submat(nfactors,0,m-1,nfactors-1).zeros();
            pt.submat(0,0,nfactors-1,nfactors-1).eye();
            at.subvec(0,nfactors-1).zeros();
            coefs.col(d+1) = at;
            PtAll.slice(d+1) = pt;
            PinfAll.slice(d+1) = pinf;
            
            //kt =     pt.cols(zind.col(0))*Z.slice(d+1).col(0);   
            //kinf = pinf.cols(zind.col(0))*Z.slice(d+1).col(0);
            for(unsigned int k = 0; k<m; k++){
              for(unsigned int l = 0; l<m1; l++){
                kt(k,0,d+1) = kt(k,0,d+1) + pt(std::max(k,zind(l,0)),std::min(k,zind(l,0)))*Z(l,0,d+1);       
              }
            }  
            for(unsigned int k = 0; k<m; k++){
              for(unsigned int l = 0; l<m1; l++){
                kinf(k,0,d+1) = kinf(k,0,d+1) + pinf(std::max(k,zind(l,0)),std::min(k,zind(l,0)))*Z(l,0,d+1);       
              }
            }  
            yhat = arma::accu(Z.slice(d+1).col(0)%at.rows(zind.col(0)));
            tmpm = kt.slice(d+1).col(0);
            zk = arma::accu(Z.slice(d+1).col(0)%tmpm.rows(zind.col(0)));     
            tmpm = kinf.slice(d+1).col(0);
            finf(0,d+1) = arma::accu(Z.slice(d+1).col(0)%tmpm.rows(zind.col(0)));
            
          }
        } else {
          
          //next observation
          //kt =  pt.cols(zind.col(j+1))*Z.slice(d).col(j+1);
          //kinf = pinf.cols(zind.col(j+1))*Z.slice(d).col(j+1);
          
          for(unsigned int k = 0; k<m; k++){
            for(unsigned int l = 0; l<m1; l++){
              kt(k,j+1,d) = kt(k,j+1,d) + pt(std::max(k,zind(l,j+1)), std::min(k,zind(l,j+1))) * Z(l,j+1,d);       
            }
          }
          for(unsigned int k = 0; k<m; k++){
            for(unsigned int l = 0; l<m1; l++){
              kinf(k,j+1,d) = kinf(k,j+1,d) + pinf(std::max(k,zind(l,j+1)), std::min(k,zind(l,j+1))) * Z(l,j+1,d);       
            }
          }
          yhat = arma::accu(Z.slice(d).col(j+1)%at.rows(zind.col(j+1)));
          tmpm = kt.slice(d).col(j+1);
          zk = arma::accu(Z.slice(d).col(j+1)%tmpm.rows(zind.col(j+1)));
          tmpm = kinf.slice(d).col(j+1);
          finf(j+1,d) = arma::accu(Z.slice(d).col(j+1)%tmpm.rows(zind.col(j+1)));
        }
        if(rankp==0){
          // d--; // negate the last increment below 
          break;
        }
      }
      
    }
    for(int i = j+1; i<p; i++){
      
      if(arma::is_finite(y(d,i))){
        ft(i,d) = zk + H(d,i);
        
        if (ft(i,d) > tiny){
          finv = 1.0/ft(i,d);
          vt(i,d) = y(d,i)-yhat;
          lik -= c + 0.5*(-log(finv) + pow(vt(i,d),2)*finv);
          at = at + kt.slice(d).col(i)*vt(i,d)*finv;
          for(unsigned int k = 0; k<m; k++){
            for(unsigned int l = k; l<m; l++){
              pt(l,k) = pt(l,k) - kt(l,i,d)*kt(k,i,d)*finv;          
            }
          }
          //pt = symmatl(pt);         
        }
      }
      if(i==(p-1)){
        //next time step         
        if(d<(n-1)){
          
          
          pt.submat(nfactors,0,m-1,nfactors-1).zeros();
          pt.submat(0,0,nfactors-1,nfactors-1).eye();
          at.subvec(0,nfactors-1).zeros();
          coefs.col(d+1)=at;
          PtAll.slice(d+1) = pt;
          //kt = pt.cols(zind.col(0))*Z.slice(d+1).col(0);   
          for(unsigned int k = 0; k<m; k++){
            for(unsigned int l = 0; l<m1; l++){
              kt(k,0,d+1) = kt(k,0,d+1) + pt(std::max(k,zind(l,0)),std::min(k,zind(l,0)))*Z(l,0,d+1);       
            }
          }   
          yhat = arma::accu(Z.slice(d+1).col(0)%at.rows(zind.col(0)));
          tmpm = kt.slice(d+1).col(0);
          zk = arma::accu(Z.slice(d+1).col(0)%tmpm.rows(zind.col(0)));         
        }
      } else {
        //next observation
        // kt = pt.cols(zind.col(i+1))*Z.slice(d).col(i+1);      
        for(unsigned int k = 0; k<m; k++){
          for(unsigned int l = 0; l<m1; l++){
            kt(k,i+1,d) = kt(k,i+1,d) + pt(std::max(k,zind(l,i+1)), std::min(k,zind(l,i+1))) * Z(l,i+1,d);       
          }
        }
        yhat = arma::accu(Z.slice(d).col(i+1)%at.rows(zind.col(i+1)));
        tmpm = kt.slice(d).col(i+1);
        zk = arma::accu(Z.slice(d).col(i+1)%tmpm.rows(zind.col(i+1))); 
        
      }
      
      
    }
    // d++;
  }    
  
  
  for(int t = d+1; t<n; t++){   
    for(int i = 0; i<p; i++){
      if(arma::is_finite(y(t,i))){
        ft(i,t) = zk + H(t,i);
        
        if (ft(i,t) > tiny){
          finv = 1.0/ft(i,t);
          vt(i,t) = y(t,i)-yhat;
          lik -= c + 0.5*(-log(finv) + pow(vt(i,t),2)*finv);
          at = at + kt.slice(t).col(i)*vt(i,t)*finv;
          for(unsigned int k = 0; k<m; k++){
            for(unsigned int l = k; l<m; l++){
              pt(l,k) = pt(l,k) - kt(l,i,t)*kt(k,i,t)*finv;          
            }
          }         
        }
      }
      
      if(i==(p-1)){
        //next time step
        if(t<(n-1)){  
          
          pt.submat(nfactors,0,m-1,nfactors-1).zeros();
          pt.submat(0,0,nfactors-1,nfactors-1).eye();
          at.subvec(0,nfactors-1).zeros();
          coefs.col(t+1)=at;
          PtAll.slice(t+1) = pt;          
          pt.submat(nfactors,0,m-1,nfactors-1).zeros();
          for(unsigned int k = 0; k<m; k++){
            for(unsigned int l = 0; l<m1; l++){
              kt(k,0,t+1) = kt(k,0,t+1) + pt(std::max(k,zind(l,0)),std::min(k,zind(l,0)))*Z(l,0,t+1);       
            }
          }        
          yhat = arma::accu(Z.slice(t+1).col(0)%at.rows(zind.col(0)));
          tmpm = kt.slice(t+1).col(0);
          zk = arma::accu(Z.slice(t+1).col(0)%tmpm.rows(zind.col(0)));         
        }
      } else {
        //next observation      
        for(unsigned int k = 0; k<m; k++){
          for(unsigned int l = 0; l<m1; l++){
            kt(k,i+1,t) = kt(k,i+1,t) + pt(std::max(k,zind(l,i+1)), std::min(k,zind(l,i+1))) * Z(l,i+1,t);       
          }
        }  
        yhat = arma::accu(Z.slice(t).col(i+1)%at.rows(zind.col(i+1)));
        
        tmpm = kt.slice(t).col(i+1);
        zk = arma::accu(Z.slice(t).col(i+1)%tmpm.rows(zind.col(i+1))); 
        
      }        
      
    } 
  }
  
  arma::mat rt(m,n,arma::fill::zeros);   
  arma::mat L0(m,m);
  //arma::mat diagm = arma::eye<arma::mat>(m,m);
  
  //This could be optimized by using symmetry!!!!
  // And computing smoothing for eta and storing only last rt0 and rt1
  // => fast state smoothing
  //and using sparsity of Z?
  
  arma::uvec ind(1);
  
  for(int t=(n-1); t>d; t--){
    ind =t;
    for(int i=(p-1); i>=0; i--){
      
      if(arma::is_finite(y(t,i)) && ft(i,t)>tiny){
        L0.eye();
        L0.cols(zind.col(i)) -= kt.slice(t).col(i)*Z.slice(t).col(i).t()/ft(i,t);
        
        rt.col(t) = L0.t()*rt.col(t);  
        rt(zind.col(i),ind) += Z.slice(t).col(i)*vt(i,t)/ft(i,t);
        
        
        
        
      }//identical(t((k[,2,drop=FALSE]%*%z[2,c(1,2,5,6),drop=FALSE]))%*%r,(t((k[,2,drop=FALSE]%*%z[2,,drop=FALSE]))%*%r)[c(1,2,5,6),,drop=FALSE])
    }
    if(t>0 && nfactors<m){
      rt.submat(nfactors,t-1,m-1,t-1)=rt.submat(nfactors,t,m-1,t);
    }
  }
  
  arma::mat rt0(m,d+1,arma::fill::zeros);
  arma::mat rt1(m,d+1,arma::fill::zeros);
  
  arma::mat Linf(m,m);   
  if(d>=0){ 
    
    
    int t=d;   
    ind =t;
    for(int i=(p-1);i>j;i--){               
      if(arma::is_finite(y(t,i)) && ft(i,t)>tiny){         
        L0.eye();
        L0.cols(zind.col(i)) -= kt.slice(t).col(i)*Z.slice(t).col(i).t()/ft(i,t);
        
        rt.col(t) = L0.t()*rt.col(t);
        rt(zind.col(i),ind) += Z.slice(t).col(i)*vt(i,t)/ft(i,t);
        
        
      }
    }
    rt0.col(d) = rt.col(d);
    
    for(int i=j;i>=0; i--){
      if(arma::is_finite(y(t,i))){
        if(finf(i,t)>tol){
          Linf.eye();
          Linf.cols(zind.col(i)) -= kinf.slice(t).col(i)*Z.slice(t).col(i).t()/finf(i,t);
          
          L0.zeros();
          L0.cols(zind.col(i)) = (ft(i,t)/finf(i,t)*kinf.slice(t).col(i)-kt.slice(t).col(i))*Z.slice(t).col(i).t()/finf(i,t);      
          
          rt1.col(t) = L0.t()*rt0.col(t) + Linf.t()*rt1.col(t); 
          rt1(zind.col(i),ind) += Z.slice(t).col(i)*vt(i,t)/finf(i,t);         
          
          rt0.col(t) = Linf.t()*rt0.col(t);
          
        } else {
          if(ft(i,t)>tiny){
            L0.eye();
            L0.cols(zind.col(i)) -= kt.slice(t).col(i)*Z.slice(t).col(i).t()/ft(i,t);        
            
            rt0.col(t) = L0.t()*rt0.col(t);
            rt(zind.col(i),ind) += Z.slice(t).col(i)*vt(i,t)/ft(i,t);
            
            rt1.col(t) = L0.t()*rt1.col(t);
            
          }
        }        
      }   
    }
    
    if(t>0 && nfactors<m){
      rt0.submat(nfactors,t-1,m-1,t-1)=rt0.submat(nfactors,t,m-1,t);
      rt1.submat(nfactors,t-1,m-1,t-1)=rt1.submat(nfactors,t,m-1,t);
    }
    
    for(int t=(d-1);t>=0;t--){
      ind=t;
      for(int i=(p-1);i>=0;i--){
        if(arma::is_finite(y(t,i))){
          if(finf(i,t)>tol){
            Linf.eye();
            Linf.cols(zind.col(i)) -= kinf.slice(t).col(i)*Z.slice(t).col(i).t()/finf(i,t);
            
            L0.zeros();
            L0.cols(zind.col(i)) = (ft(i,t)/finf(i,t)*kinf.slice(t).col(i)-kt.slice(t).col(i))*Z.slice(t).col(i).t()/finf(i,t);      
            
            rt1.col(t) = L0.t()*rt0.col(t) + Linf.t()*rt1.col(t); 
            rt1(zind.col(i),ind) += Z.slice(t).col(i)*vt(i,t)/finf(i,t);         
            
            rt0.col(t) = Linf.t()*rt0.col(t);            
            
            
          } else {
            if(ft(i,t)>tiny){
              L0.eye();
              L0.cols(zind.col(i)) -= kt.slice(t).col(i)*Z.slice(t).col(i).t()/ft(i,t);        
              
              rt0.col(t) = L0.t()*rt0.col(t);
              rt(zind.col(i),ind) += Z.slice(t).col(i)*vt(i,t)/ft(i,t);
              
              rt1.col(t) = L0.t()*rt1.col(t);
            }
          }        
        }   
      }
      if(t>0 && nfactors<m){
        rt0.submat(nfactors,t-1,m-1,t-1)=rt0.submat(nfactors,t,m-1,t);
        rt1.submat(nfactors,t-1,m-1,t-1)=rt1.submat(nfactors,t,m-1,t);
      }
    }
  }
  //Lazy solution (uses symmatl):
  
//    Rcpp::Rcout<<"coefs.row(2) "<<coefs.row(2)<<std::endl;
  for(int t=0; t<=d; t++){     
    PtAll.slice(t) = symmatl(PtAll.slice(t));
    PinfAll.slice(t) = symmatl(PinfAll.slice(t));
    
    coefs.col(t) = coefs.col(t) + PtAll.slice(t)*rt0.col(t) + PinfAll.slice(t)*rt1.col(t);    
  }
  for(int t=d+1; t<n;t++){
    PtAll.slice(t) = symmatl(PtAll.slice(t));
    coefs.col(t) = coefs.col(t) + PtAll.slice(t)*rt.col(t);
  }
  for(int t = 0; t<n; t++){
    ind =t;
    for(int i = 0; i<p; i++){
      theta(t,i) = arma::accu(Z.slice(t).col(i)%coefs(zind.col(i),ind));
    }
  }  
  
  return lik;
}   


