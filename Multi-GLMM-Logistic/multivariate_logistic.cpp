#include <TMB.hpp>

template<class Type>
Type objective_function<Type>::operator() ()
{
  using namespace density;
  
  DATA_MATRIX(X);             // n x p
  DATA_MATRIX(Y);             // n x r
  DATA_MATRIX(Ntrials);       // n x r
  DATA_IVECTOR(hospital_id);  // H
  
  PARAMETER_MATRIX(beta);     // p x r  
  PARAMETER_MATRIX(U);        // H x r
  PARAMETER_VECTOR(rho);      // r(r-1)/2
  PARAMETER_VECTOR(sigma);    // r
  
  parallel_accumulator<Type> nll(this);
  
  const int n = Y.rows();
  const int p = X.cols();
  const int r = Y.cols();
  const int H = U.rows(); 
  
  matrix<Type> Xbeta = X * beta;  // n x r
  
  // likelihood: binomial with logit link
  for (int j = 0; j < r; j++) {
    for (int i = 0; i < n; i++) {
      int h = hospital_id(i);     
      Type eta = Xbeta(i, j) + U(h, j);
      Type pj  = invlogit(eta);     
      nll -= dbinom(Y(i, j), Ntrials(i, j), pj, true);
    }
  }
  
  // random effects: U ~ MVN(0, Sigma)
  UNSTRUCTURED_CORR_t<Type> Corr(rho);
  for (int h = 0; h < H; h++) {
    nll += VECSCALE(Corr, sigma)(U.row(h));
  }
  
  matrix<Type> Cor = Corr.cov();        
  ADREPORT(Cor);
  vector<Type> SigmaDiag = sigma;     
  ADREPORT(SigmaDiag);
  
  return nll;
}
