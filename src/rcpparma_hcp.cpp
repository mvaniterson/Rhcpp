// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List rcpparma_hcp(NumericMatrix Fr, NumericMatrix Yr, int k, int lambda1, int lambda2, int lambda3, int iter)
{
  double tol = 1e-6;
  int ii = 0;
  int n1 = Fr.nrow();
  int d1 = Fr.ncol();
  int n2 = Yr.nrow();
  int d2 = Yr.ncol();

  arma::mat F(Fr.begin(), n1, d1, false);
  arma::mat Y(Yr.begin(), n2, d2, false);

  /*int k = Rcpp::as<int>(kr);
  int iter = Rcpp::as<int>(iterr);
  int lambda1 = Rcpp::as<int>(lambda1r);
  int lambda2 = Rcpp::as<int>(lambda2r);
  int lambda3 = Rcpp::as<int>(lambda3r);*/

  Rcpp::NumericMatrix Ur(d1, k);
  Rcpp::NumericMatrix Zr(n1, k);
  Rcpp::NumericMatrix Br(k, d2);

  Rcpp::NumericVector unif = runif(k*d2);
  int cntr = 0;
  for(int i=0; i<k; i++)
    for(int j=0; j<d2; j++)
      {
        Br[i,j] = unif[cntr];
        cntr += 1;
      }

  Rcpp::NumericVector o(iter);

  arma::mat U = Rcpp::as<arma::mat>(Ur);
  arma::mat Z = Rcpp::as<arma::mat>(Zr);
  arma::mat B = Rcpp::as<arma::mat>(Br);

  arma::mat diagB = arma::eye<arma::mat>(k, k);
  arma::mat diagZ = arma::eye<arma::mat>(k, k);
  arma::mat diagU = arma::eye<arma::mat>(d1, d1);

  if(iter > 0)
    {
      for (ii=0; ii < iter; ii++) {
        /*
          o[ii] = sqrt(trace((Y-Z*B)*(Y-Z*B).t())) + sqrt(trace((Z-F*U)*(Z-F*U).t()))*lambda1 + lambda2*sqrt(trace(B*B.t())) + lambda3*sqrt(trace(U*U.t()));
          Z = Y*trans(B) + lambda1*F*U*inv(B*trans(B) + lambda1*diagB);
          B = solve(Z.t()*Z + lambda2*diagZ, Z.t()*Y);
          U = solve(F.t()*F*lambda1 + lambda3*diagU, lambda1*F.t()*Z);
        */

        o[ii] = norm(Y-Z*B,2) + norm(Z-F*U,2)*lambda1 + norm(B,2)*lambda2 + lambda3*norm(U,2);
        //o[ii] = accu(pow(Y-Z*B,2)) + accu(pow(Z-F*U,2))*lambda1 + accu(pow(B,2))*lambda2 + lambda3*accu(pow(U,2));
        Z = (Y*B.t() + lambda1*F*U)*(B*B.t() + lambda1*diagB).i();
        B = (Z.t()*Z + lambda2*diagZ).i()*(Z.t()*Y);
        U = (F.t()*F + (lambda3/lambda1)*diagU).i()*(F.t()*Z);

        if(ii > 0) {
          if((abs(o[ii] - o[ii-1])/o[ii]) < tol)
            break;
        }
      }
    }
  
  return Rcpp::List::create(Rcpp::Named("Z")=Z,Rcpp::Named("B")=B,Rcpp::Named("U")=U,Rcpp::Named("o")=o, Rcpp::Named("iter")=ii);
}

