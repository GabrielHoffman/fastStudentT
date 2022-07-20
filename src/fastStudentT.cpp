// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"


// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// Fast computation of Student t log likelihood function that removes redundant computation
//
// [[Rcpp::export]]
double dt_fast(const arma::vec &y, double mu, double s, double nu){
    
    int n = y.size();
    arma::vec x = (y-mu) / s;
    
    double value = n*lgamma((nu+1)/2); 
    value -= n*0.5*log(nu*M_PI); 
    value -= n*lgamma(nu/2);
    arma::vec tmp = pow(x,2);
    value -= (nu+1)/2 * sum(log1p(tmp / nu));

    return value;
}

// Fast estimation of mu and s for Student t, given df
//
// [[Rcpp::export]]
Rcpp::List studentTEM( const arma::vec &y, double nu, int niter){
  
    int n = y.size();
    arma::vec w(n);
    arma::vec tmp; 
    w.fill(1.0);
    double mu, sigSq;

    double mu_prev = 0, sigSq_prev = 0;
    int j;

    // Estimate mu and sigSq
    for(j=0; j<niter; j++){
        // Fit EM algorithm
        mu = as_scalar(y.t()*w) / as_scalar(w.t()*w);
        tmp = pow(y - mu,2);
        sigSq = as_scalar(w.t()*tmp) / (double) n;
        w = ((nu + 1) * sigSq) / (nu*sigSq + tmp);

        // check convergence
        if( (fabs(mu_prev - mu) < 1e-7) && (fabs(sigSq - sigSq_prev) < 1e-7) ){
            break;
        } 

        // store values for convergence check   
        mu_prev = mu;
        sigSq_prev = sigSq;
    }

    return Rcpp::List::create(Rcpp::Named("m") = mu, 
                        Rcpp::Named("s") = sqrt(sigSq),
                        Rcpp::Named("iter") = j);
}