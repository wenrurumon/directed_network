
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <unistd.h>

//#include "ADMM_methods.cpp"
#include <math.h>

using namespace arma;

#include <math.h>
using namespace arma;



inline double covsel_obj(const mat & s, const mat & x, const mat & z, const double & lambda){
  return( trace(s*x) -log(det(x)) + lambda*accu(abs(z)) );
}


double norm_fro(const mat & x){
  return( sqrt(trace( (x.t() * x) )) );
}


// [[Rcpp::export]]
mat shrinkage(const mat & a, double kappa) {
  /*
  Rcpp::Rcout << "kappa " << kappa << std::endl
              << "a.max " << a.max() << std::endl
              << "a.min " << a.min() << std::endl;
  */
  mat s1 = clamp(  a-kappa, 0,  std::max( a.max()-kappa, kappa));
  mat s2 = clamp( -a-kappa, 0,  std::max(-a.min()-kappa, kappa));
  s1 = s1 - s2;
  return( s1 );
}

// [[Rcpp::export]]
mat covsel(const mat & d, const int max_iter=1000, double alpha=1.2,
           double lambda=0.0001, double rho=0.5){
  double ABSTOL = 1e-5;
  double RELTOL = 1e-3;
  mat s = cov(d);
  unsigned int n = s.n_rows;
  mat x = zeros(n,n);
  mat z = zeros(n,n);
  mat u = zeros(n,n);

  double  r_norm = 0, s_norm = 0,  eps_pri = 0, eps_dual = 0;
  //double objval;
  for( int i = 0 ; i < max_iter && r_norm <= eps_pri && s_norm <= eps_dual; i++){
    vec es;
    mat Q;
    eig_sym(es,Q,rho*(z-u) - s);
    vec xi = (es +sqrt(pow(es,2)+4*rho))/(2*rho);
    x = Q * diagmat(xi) * Q.t();

    mat zold = z;
    mat x_hat = alpha*x + (1-alpha)*zold;
    z = shrinkage( x_hat +u, lambda/rho);

    u = u + (x_hat - z);

    //objval =  covsel_obj(s,x,z,lambda);
    r_norm = norm_fro( x - z);
    s_norm = -rho*norm_fro(z - zold);

    eps_pri = n*ABSTOL + RELTOL*std::max(norm_fro(x), norm_fro(z) );
    eps_dual = n*ABSTOL + rho*RELTOL*norm_fro(u);
    /*
    Rcpp::Rcout << "i " << i  << std::endl;
    Rcpp::Rcout << "r_norm " << r_norm <<std::endl;
    Rcpp::Rcout << "s_norm " << s_norm <<std::endl;
    Rcpp::Rcout << "eps_pri " << eps_pri <<std::endl;
    Rcpp::Rcout << "eps_dual " << eps_dual <<std::endl;
    Rcpp::Rcout << "objval " << objval <<std::endl;
    */
  }
  return x;

}

void factor(const mat & A, mat & L, const double & rho){
  unsigned m = A.n_rows;
  unsigned n = A.n_cols;

  if ( m >=n ){
    chol(L, A.t()*A+rho*eye(n,n), "lower");
  }else{
    chol(L, eye(m,m)+1/rho*(A*A.t()),"lower");
  }
}

double norm_1(const mat & x){
  return( as_scalar( max(sum(abs(x) ) ) ) );
}

double norm(const mat & x){
  return( as_scalar( max(svd(x) ) ) );
}

double objective_lasso(const mat & A, const vec & b, const vec & x,
                       const vec & z, double lambda=0.0001){
  return( 0.5*sum(pow(A*x - b, 2))+ lambda*norm_1(z) );
}


// [[Rcpp::export]]
vec lasso(const mat & A, const vec & b, const int max_iter=1000, double alpha=1.2,
          double lambda=0.1, double rho=1){
  double ABSTOL = 1e-4;
  double RELTOL = 1e-2;
  unsigned int m = A.n_rows;
  unsigned int n = A.n_cols;

  mat At = A.t();
  mat Atb = At*b;
  mat L;

  factor(A,L,rho);
  mat U = L.t();

  vec x = zeros(n);
  vec z = zeros(n);
  vec u = zeros(n);

  vec q;
  double  r_norm = 0, s_norm = 0, eps_pri = 0, eps_dual = 0;
  //double objval;
  double sqrt_n = sqrt(n);
  double rho_2 = rho*rho;
  for( int i = 0 ; i < max_iter && r_norm <= eps_pri && s_norm <= eps_dual; i++){
    q = Atb + rho *(z-u);

    if( m >= n ){
      solve(x,U,solve(L,q));
    }else{
      x = q/rho-(At*(solve(U,solve(L, A*q))))/rho_2;
    }

    vec zold = z;
    vec x_hat = alpha*x + (1-alpha)*zold;
    z = shrinkage(x_hat+u,lambda/rho);

    u = u+(x_hat -z);
    //objval = objective_lasso(A,b,x,z,lambda);

    r_norm = norm( x - z);
    s_norm = -rho*norm(z - zold);

    eps_pri = sqrt_n*ABSTOL + RELTOL*std::max(norm(x), norm(-z) );
    eps_dual = sqrt_n*ABSTOL + rho*RELTOL*norm(u);
    /*
     Rcpp::Rcout << "i " << i  << std::endl;
     Rcpp::Rcout << "r_norm " << r_norm <<std::endl;
     Rcpp::Rcout << "s_norm " << s_norm <<std::endl;
     Rcpp::Rcout << "eps_pri " << eps_pri <<std::endl;
     Rcpp::Rcout << "eps_dual " << eps_dual <<std::endl;
     Rcpp::Rcout << "objval " << objval <<std::endl;
     */
  }
  return z;

}




// [[Rcpp::export]]
arma::mat ca_mean(const arma::mat & x){
  arma::mat m = sum(x)/x.n_rows;
  return(m);
}


// [[Rcpp::export]]
arma::mat ca_var(const arma::mat & x){
  return(var(x));
}



arma::mat center(const arma::mat & x){
  return( x - repmat(sum(x)/x.n_rows,x.n_rows,1) );
}

void center_replace( arma::mat & x){
  x = x - repmat(sum(x)/x.n_rows,x.n_rows,1);
}

void scale_replace( arma::mat & x){
  x = x / repmat(sqrt(var(x)),x.n_rows,1);
}

arma::mat scale( arma::mat & x){
  return( x / repmat(sqrt(var(x)),x.n_rows,1) );
}

// [[Rcpp::export]]
arma::mat center_scale(const arma::mat & x){
  arma::mat m = sum(x)/x.n_rows;
  return( (x -repmat(m,x.n_rows,1)) / repmat(sqrt(var(x)),x.n_rows,1) ); ////TODO: add check if any var = 0
}

void center_scale_replace(arma::mat & x){
  arma::mat m = sum(x)/x.n_rows;
  x = (x -repmat(m,x.n_rows,1)) / repmat(sqrt(var(x)),x.n_rows,1) ; ////TODO: add check if any var = 0
}



// [[Rcpp::export]]
mat initial_sem_nofix( const mat & data, const uvec & variable,
                double rho = 1.0, double lambda = 0.1, double prec=1e-5, int max_iter = 1000){
  unsigned int n = data.n_rows;
  unsigned int d = data.n_cols - 1;

  int num_scr;
  if( d>n){
    if( n <=3 ){
      num_scr = n;
    }else{
      num_scr = std::ceil( n/log(n) );
    }
  }else{
    if(d<=4){
      num_scr = d;
    }else{
      num_scr = std::ceil( sqrt(d) );
    }
  }

  mat td = center_scale(data) ;
  mat coef_mat( td.n_cols, td.n_cols );
  coef_mat.fill(0.0);


  for( unsigned int index = 0; index < variable.n_rows; index++){
    unsigned int idx = variable(index) - 1;
    mat W = td;
    W.shed_col( idx  );
    mat Y = td.col(idx );
    vec xy0 = W.t() * Y;
    uvec order = sort_index(xy0,"descend");
    uvec idx_scr = order.subvec(0,num_scr-1);
    W = W.cols(idx_scr);
    idx_scr.elem( find( idx_scr>=idx )  ) += 1;

    vec beta = lasso(W,Y);
    //Rcpp::Rcout << "beta " << beta << std::endl;
    for(unsigned int i = 0; i < idx_scr.n_rows; i++ ){
      coef_mat(idx_scr(i),idx) = beta(i);
    }
    //Rcpp::Rcout << "idx_scr " << idx_scr << std::endl;
  }

  return( coef_mat.cols(variable-1).t() );
}


// [[Rcpp::export]]
mat initial_sem(const mat & data, const uvec &  fixed, const uvec & variable,
                double rho = 1.0, double lambda = 0.1, double prec=1e-5, int max_iter = 1000){
  mat x = data.cols( fixed - 1);
  mat xt = x.t();
  mat y = data.cols(variable - 1);

  mat y_hat = x* covsel(xt*x)*xt*y;

  return( initial_sem_nofix( join_rows(y_hat, x),  variable) );
}

mat initial_sem_stable(const uvec &  phenos, const uvec &  fixed, const mat & data,
                double rho = 1, int max_iter = 100, double lambda=0.1){
  unsigned int n = data.n_rows;
  unsigned int d = data.n_cols - 1;
  //unsigned int maxdf = (n<d)?d:n;
  mat data_var = ca_var(data);
  mat csx = center_scale(data);

  int num_scr1, num_scr2;
  if( d>n){
    if( n <=3 ){
      num_scr1 = n;
      num_scr2 = n;
    }else{
      num_scr1 = std::ceil( n/log(n) );
      num_scr2 = n-1;
    }
  }else{
    if(d<=4){
      num_scr1 = d;
      num_scr2 = d;
    }else{
      num_scr1 = std::ceil( sqrt(d) );
      num_scr2 = std::ceil( d/log(d) );
    }
  }

  for( unsigned int index = 0; index < phenos.n_rows; index++){
    mat x0 = csx;
    x0.shed_col(phenos(index) - 1);
    mat y = csx.unsafe_col(phenos(index) - 1);
    vec xy0 = x0.t() * y;
    mat beta(d-1,1);
    uvec order0 = sort_index(xy0,"descend");
    uvec idx_scr = order0;
    //unsigned int num_scr = idx_scr.n_rows;
    uvec idx_scr1 = order0.subvec(0,num_scr1-1);
    uvec idx_scr2 = order0.subvec(0,num_scr2-1);
    mat x1 = x0.cols(idx_scr);
    mat xx = x1.t() * x1;
    //double gamma = sum(abs(xx)).max();
    vec beta0;
    beta0.fill(0);


  }
  return(csx);
}

void lasso_ladm_scr(const mat & xy0, const mat & x0, const mat & xx0, double & T,
                    mat & xy, mat & x, mat & xx, double rho, int max_ite,
                    double prec, const vec & beta, const uvec & idx_scr,
                    double lambda, bool flag ){
  vec beta0, beta1;
  unsigned int n = idx_scr.n_cols;
  beta0 = beta.rows(idx_scr);


  uvec idx_a(n), idx_i(n), idx_a1(n), idx_i1(n);
  idx_i.fill(1);
  idx_i.elem( find(beta0!=0)).zeros();

  for( unsigned int j = 0,size_a = 0; j < n; j++){
    if ( beta0[j] != 0){
      idx_a[size_a] = j;
      size_a++;
    }
  }

  x = x0.cols(idx_scr);
  xy = xy0.cols(idx_scr);
  if(flag) {
    xx = xx0;
    T = accu(abs(xx))*0.8;
  }
}

void lineaization_lasso(vec & beta0, vec & beta1, vec & beta_tild, vec & g,
                        uvec & idx_a1, uvec & idx_i1, double T,
                        double threshold, int dim){
  int size_a1 = 0;
  for( int j = 0; j < dim; j++){
    beta_tild(j) = beta0(j) - g(j)/T;
    if( fabs(beta_tild[j]) <= threshold ){
      beta1(j) = 0;
      idx_i1(j) = 1;
    }else{
      idx_a1(size_a1) = j;
      size_a1++;
      idx_i1(j) = 0;
      if( beta_tild(j) > 0 ){
        beta1(j) = beta_tild(j) - threshold;
      }else{
        beta1(j) = beta_tild(j) + threshold;
      }
    }
  }
}
