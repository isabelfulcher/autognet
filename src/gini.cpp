#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// Sort vector
NumericVector stl_sort(NumericVector x) {
   NumericVector y = clone(x);
   std::sort(y.begin(), y.end());
   return y;
}

// Remove NAs from vector
NumericVector naomit(NumericVector x){
  std::vector<double> r(x.size());
  int k=0;
    for (int i = 0; i < x.size(); ++i) {
      if (x[i]==x[i]) {
        r[k] = x[i];
        k++;
      }
    }
 r.resize(k);
 return wrap(r);
}

// Compute weighted ( by index )  sum of ordered vector
double osum (NumericVector x){
  int n = x.length();
  double total = 0;
  for (int i = 0; i < n ; ++i) {
    total += ((i+1)*x(i));
  }
  return total;
}

//' Calculate the Gini Index using Rcpp
//'
//' Given a vector of numbers, calcuate the
//' Gini index being sensitive to NAs
//'
//' @param x A numeric vector
//' @return A float of the Gini index where values are
//' between [0,1]
//'
//' @examples
//'
//' x <- runif(1000)
//'
//' @export
// [[Rcpp::export]]
float giniCpp (NumericVector x){

  int n = x.length();
  NumericVector sortedrow = stl_sort(naomit(x));
  float G = ((2*osum(sortedrow))/sum(sortedrow)) - (n + 1);
  float gini = G/n;
  return (gini);

}
