// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List rcpparma_hcpp(NumericMatrix Fr, NumericMatrix Yr, NumericMatrix xr, int k, int lambda1, int lambda2, int lambda3, int iter)
{
  double tol = 1e-6;
  double xtx;
  int ii = 0;
  int n1 = Fr.nrow();
  int d1 = Fr.ncol();
  int n2 = Yr.nrow();
  int d2 = Yr.ncol();

  arma::mat F(Fr.begin(), n1, d1, false);
  arma::mat Y(Yr.begin(), n2, d2, false);
  arma::mat x(xr.begin(), n1, 1, false);

  /*int k = Rcpp::as<int>(kr);
  int iter = Rcpp::as<int>(iterr);
  int lambda1 = Rcpp::as<int>(lambda1r);
  int lambda2 = Rcpp::as<int>(lambda2r);
  int lambda3 = Rcpp::as<int>(lambda3r);*/

  Rcpp::NumericMatrix Ur(d1, k);
  Rcpp::NumericMatrix Zr(n1, k);
  Rcpp::NumericMatrix Br(k, d2);
  Rcpp::NumericMatrix gr(1, d2);

  Rcpp::NumericVector unif = runif(k*d2);
  //int cntr = 0;
  for(int i=0; i<k; i++) 
    for(int j=0; j<d2; j++)
      {
        //Br[i,j] = unif[cntr];
        //cntr += 1;
        Br[i,j] = (i+j)/(k*d2);
      }

    for(int j=0; j<d2; j++)
      gr[1,j] = 0;

  Rcpp::NumericVector o(iter);

  arma::mat U = Rcpp::as<arma::mat>(Ur);
  arma::mat Z = Rcpp::as<arma::mat>(Zr);
  arma::mat B = Rcpp::as<arma::mat>(Br);
  arma::mat g = Rcpp::as<arma::mat>(gr);

  arma::mat diagB = arma::eye<arma::mat>(k, k);
  arma::mat diagZ = arma::eye<arma::mat>(k, k);
  arma::mat diagU = arma::eye<arma::mat>(d1, d1);

  xtx = arma::as_scalar(x.t()*x);

  if(iter > 0)
    {
      for (ii=0; ii < iter; ii++) {

        o[ii] = pow(norm(Y-x*g-Z*B),2) + pow(norm(Z-F*U),2)*lambda1 + pow(norm(B),2)*lambda2 + lambda3*pow(norm(U),2); //for matrices Frobenius norm is default p = 2 or p = "fro"

        Z = ((Y - x*g)*B.t() + lambda1*F*U)*(B*B.t() + lambda1*diagB).i();

        B = solve(Z.t()*Z + lambda2*diagZ, Z.t()*(Y-x*g));

        U = solve(F.t()*F + (lambda3/lambda1)*diagU, F.t()*Z);

        g = x.t()*(Y-Z*B)/xtx;

        if(ii > 0) {
          if((abs(o[ii] - o[ii-1])/o[ii]) < tol)
            break;
        }
      }
    }
  
  return Rcpp::List::create(Rcpp::Named("Z")=Z,Rcpp::Named("B")=B,Rcpp::Named("U")=U,Rcpp::Named("g")=g,Rcpp::Named("o")=o,Rcpp::Named("iter")=ii);
}

