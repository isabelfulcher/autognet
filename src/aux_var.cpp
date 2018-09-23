#include <math.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using namespace Rcpp;

int rcpp_rbinom_one(float prob) {
  int v = as<int>(rbinom(1, 1, prob));
  return(v);
}


// [[Rcpp::export]]
IntegerVector oneMultinomCall(NumericVector probs) {
  int k = probs.size();
  IntegerVector ans(k);
  rmultinom(1, probs.begin(), k, ans.begin());
  return(ans);
}

// [[Rcpp::export]]
IntegerVector oneMultinomC(NumericVector probs) {
  int k = probs.size();
  SEXP ans;
  PROTECT(ans = Rf_allocVector(INTSXP, k));
  probs = Rf_coerceVector(probs, REALSXP);
  rmultinom(1, REAL(probs), k, &INTEGER(ans)[0]);
  UNPROTECT(1);
  return(ans);
}

// [[Rcpp::export]]
IntegerVector callRMultinom(NumericVector x) {
  int n = x.size();
  IntegerVector d(n);
  R::rmultinom(1, x.begin(), n, d.begin());
  return d;
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
          // Rcpp::Rcout << i;
          // Multinomial case
          NumericVector prob_vec(group_length + 1);
          prob_vec = prob_vec + 1;
          for (int m = 0; m < group_length; ++m){
            int jprime = j + m; //already set j to 0 above whereas it is 1 in R
            arma::vec rowVec = rho_mat(_,j);
            arma::mat covVec = vectorise(cov_mat.rows(i,i));
            arma::uvec whichN = adjacency[i];
            arma::mat covAdjVec = cov_mat.cols(j,j);
            float weights_i = weights[i];
            float nei_w = sum(covAdjVec.elem(whichN)/weights_i);
            float prob_Lj = exp(arma::as_scalar(tau[jprime] + dot(rowVec,covVec) + nu[jprime]*nei_w));
            prob_vec(m) = prob_Lj;
          }
          // make the rmultinom call
          NumericVector prob_vec2 = prob_vec / sum(prob_vec);
          IntegerVector multi_out = callRMultinom(prob_vec2);

          // update elements in matrix via another loop
          for (int m = 0; m < group_length; ++m){
            // float x = j+m;
            // Rcpp::Rcout << x;
            cov_mat(i,j + m) = multi_out(m);
          }
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

        } // add in normal here once form is decided

        j = j + group_length;

      }
    }
  }
  return(cov_mat);

}

