#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::depends("RcppArmadillo")]]

// [[Rcpp::export]]
double score_function_regression_with_fix(int node, const uvec & parents, const uvec & fixed, const mat & data){
  mat W = data.cols( join_vert(parents-1,fixed-1) );
  mat X = data.cols( fixed );
  mat Y = data.unsafe_col(node-1);
  mat xtx = X.t() * X;
  mat xtx_inv = pinv( (xtx.t() + xtx )/2.0 );
  mat wtx = W.t() * X;
  mat xtw = wtx.t();
  mat wtxxtx_inv_xt =  wtx * xtx_inv * X.t();
  mat tmp = wtxxtx_inv_xt * W;
  mat delta = pinv( (tmp+tmp.t())/2 ) * wtxxtx_inv_xt * Y;
  
  mat residual = pow(Y - W*delta,2);
  return( accu(residual) );
}

// [[Rcpp::export]]
double score_function_regression(int node, const uvec & parents, const mat & data){
  mat X = data.cols( parents - 1 );
  mat Y = data.unsafe_col(node - 1);
  vec coef = solve(X,Y);  
  mat residual = pow(Y - X*coef,2);
  return( accu(residual) );
}