#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// function of posterior mean
double mu_post(const arma::colvec x, double x_sigma, double mu_prior, double sigma_prior){
  int n = x.n_rows;
  double a = 1/(1/pow(sigma_prior,2) + n/pow(x_sigma,2));
  double b = mu_prior/pow(sigma_prior,2) + sum(x)/pow(x_sigma,2);
  return a*b;
  }

// function of posterior variance
double sigma_post(const arma::colvec x, double x_sigma, double sigma_prior){
  int n = x.n_rows;
  double a = 1/pow(sigma_prior,2) + n/pow(x_sigma, 2);
  return 1/a;
}

// [[Rcpp::export]]
arma::mat post_rsmp(const arma::colvec yn, double y_sigma, double mu_prior, double sigma_prior, const int B, const int N) {
  arma::colvec y;
  arma::mat yN(B, N);
  arma::mat z(B, N, fill::randn);
  arma::mat mu_seq(B, N+1);
  arma::mat sigma_seq(B, N+1);
  
  for(int j = 0; j < B; j++){
    y = yn;
    mu_seq(j,0) = mu_post(yn, y_sigma, mu_prior, sigma_prior);
    sigma_seq(j,0) = sqrt(sigma_post(yn, y_sigma, sigma_prior) + pow(y_sigma, 2));
    
    for(int i = 0; i < N; i++){
      yN = mu_seq(j,i) + z(j,i) * sigma_seq(j,i);
      y.insert_rows(y.n_rows, yN);
      mu_seq(j,i+1) = mu_post(y, y_sigma, mu_prior, sigma_prior);
      sigma_seq(j,i+1) = sqrt(sigma_post(y, y_sigma, sigma_prior) + pow(y_sigma, 2));
      }
  }
  return mu_seq;
}


