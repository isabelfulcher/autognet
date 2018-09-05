#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

int rcpp_rbinom_one(float prob) {
  int v = as<int>(rbinom(1, 1, prob));
  return(v);
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
NumericMatrix auxVarCpp (NumericVector tau, NumericVector rho, NumericVector nu,
                         int N, int R, int J, NumericMatrix rho_mat,
                         List adjacency, NumericMatrix cov_i, IntegerVector weights,
                         IntegerVector group_lengths, IntegerVector group_functions){

  NumericMatrix cov_mat = cov_i;
  int number_of_groups = group_lengths.size();

  // Number of iterations
  for (int r = 0; r < R; ++r){

    // Number of people
    for (int i = 0; i < N; ++i){

      // Index covariate
      int j = 0;

      // Number of groups (of covariants)
      for (int group_index = 0; group_index < number_of_groups; ++group_index){

        // Values associated with each group
        int group_length = group_lengths[group_index];
        int group_function = group_functions[group_index];

        // group_length is the number of binarized covariates
        // group_function is a numeric dictionary key that utilizes a value from the covariate_process function
        // if group_length is > 1, meaning that we have multiple binarized covariates associated with a specific
        // entity, then it is definitely a multinomial

        if(group_length > 1){

          // Multinomial case
          // prob_vec <- sapply(0:(group_length -1), function(m){
          //   j_prime <- j + m
          //   exp(tau[j_prime] + sum(rho_mat[,j_prime]*cov.mat[i,]) + nu[j_prime]*sum(cov.mat[adjacency[[i]],j_prime]/weights[i]))
          // })

          // for the rmultinom call, have to append a 1 and remove the last value; update several values
          // cov.mat[i, j + (0:(group_length -1))] <- (rmultinom(1,1,c(prob_vec,1))[,1])[-1*group_length]

        } else if(group_function == 1){

          // Logistic / binary case

          // Fix this probably have to use a loop because vectorized logic is hard here
          NumericVector rowVec = rho_mat(_,j);
          NumericVector covVec = cov_mat(i,_);
          NumericVector whichN = adjacency[i];
          NumericVector covAdjVec = cov_mat(whichN,j)
          float prob_Lj = plogis(tau[j] + sum(rowVec*covVec)) + nu[j]*sum(covAdjVec/weights[i]));
          cov_mat(i,j) = rcpp_rbinom_one(0.5); // prob_Lj

        } // add in normal here once form is decide

        j = j + group_length;

      }
    }
  }
  return(cov_mat);

}

