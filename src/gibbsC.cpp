#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs sampler using Rcpp
//' @description A Gibbs sampler using Rcpp
//' @param N the length of the chain
//' @param thin the number of Binomial(thin, y)
//' @return a random sample of size \code{n}
//' @examples
//' \dontrun{
//' rnC <- gibbsC(1000,50)
//' plot(rnC[1,],type='l')
//' plot(rnC[2,],type='l')
//' }
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int thin) {
  NumericMatrix X(2,N);
  NumericVector temp(1);
  X(0,0) = 0;
  X(1,0) = 0;
  double a = 1, b =2;
  double X1 = 0;
  double Y1 =0;
  for(int i=1;i<=N;i++){
    Y1 =  X(1,i-1);
    temp =  rbinom(1,thin,Y1);
    X(0,i) = temp[0];
    X1 = X(0,i);
    temp = rbeta(1,X1+a,thin-X1+b);
    X(1,i) = temp[0];
  }
  return(X);
}
