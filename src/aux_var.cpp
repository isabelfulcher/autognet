#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

int rcpp_rbinom_one(float prob) {
  int v = as<int>(rbinom(1, 1, prob));
  return(v);
}


//' Run Gibbs sampler for auxiliary covariate values using Rcpp
//'
//' Given the specific inputs, determine auxiliary covariate
//' values using a Gibbs sampling procedure.
//'
//' @param tau A numeric vector for the intercept terms in the covariate model
//' @param rho A numeric vector for the correlation terms in the covariate model
//' @param nu A numberic vector for the neighbor terms in the covariate model
//' @param N An integer indicating the size of the interconnected network
//' @param R An integer indicating the number of iterations for the Gibbs
//' @param J An integer for the number of covariates
//' @param rho_mat A numeric matrix for rho terms
//' @param adjacency A binary matrix indicating connected units
//' @param cov_i A numeric matrix for observed covariate values (starting values for chain)
//' @param weights A numeric vector indicating the number of neighbors for each node
//' @param group_lengths An integer vector indicating the number of categories for each variable
//' @param group_functions An integer vector indicating the type of variable
//' @return A numeric matrix for auxiliary covariate values
//' between [0,1]
//'
//' @examples
//'
//' x <- runif(1000)
//'
//' @export
// [[Rcpp::export]]
arma::mat auxVarCpp (NumericVector tau, NumericVector rho, NumericVector nu,
                     int N, int R, int J, NumericMatrix rho_mat,
                     List adjacency, arma::mat cov_i, IntegerVector weights,
                     IntegerVector group_lengths, IntegerVector group_functions){

  arma::mat cov_mat = cov_i;
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
          arma::vec rowVec = rho_mat(_,j);
          arma::mat covVec = vectorise(cov_mat.rows(i,i));
          arma::uvec whichN = adjacency[i];
          arma::mat covAdjVec = cov_mat.cols(j,j);
          float weights_i = weights[i];
          float nei_w = sum(covAdjVec.elem(whichN)/weights_i);
          float prob_Lj = R::plogis(arma::as_scalar(tau[j] + dot(rowVec,covVec) + nu[j]*nei_w), 0, 1, 1, 0);
          cov_mat(i,j) = rcpp_rbinom_one(prob_Lj); // prob_Lj

        } // add in normal here once form is decide

        j = j + group_length;

      }
    }
  }
  return(cov_mat);

}

