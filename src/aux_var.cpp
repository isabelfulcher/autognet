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

          // Multinomial case
          NumericVector prob_vec(group_length + 1, 1.0);
          for (int m = 0; m < group_length; ++m){
            int j_prime = j + m; //already set j to 0 above whereas it is 1 in R
            arma::vec rowVec = rho_mat(_,j_prime);
            arma::mat covVec = vectorise(cov_mat.rows(i,i));
            arma::uvec whichN = adjacency[i];
            arma::mat covAdjVec = cov_mat.cols(j_prime,j_prime);
            float weights_i = weights[i];
            float nei_w = sum(covAdjVec.elem(whichN)/weights_i);
            float prob_Lj = exp(arma::as_scalar(tau[j_prime] + dot(rowVec,covVec) + nu[j_prime]*nei_w));
            prob_vec(m) = prob_Lj;
          }

          // make the rmultinom call
          NumericVector prob_vec2 = prob_vec / sum(prob_vec);
          IntegerVector multi_out = callRMultinom(prob_vec2);
          //Rcpp::Rcout << prob_vec2 << "\n";
          // update elements in matrix via another loop
          for (int m = 0; m < group_length; ++m){
            cov_mat(i,j + m) = multi_out(m);
          }
        } else if(group_function == 1){

          // Logistic / binary case

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

//' Run Gibbs sampler for auxiliary outcome values using Rcpp
//'
//' Given the specific inputs, determine auxiliary outcome
//' values using a Gibbs sampling procedure.
//'
//' @param beta A numeric vector of parameters from outcome model
//' @param trt A numeric vector of the treatment values
//' @param cov A numeric matrix for observed covariate values (starting values for chain)
//' @param N An integer indicating the size of the interconnected network
//' @param R An integer indicating the number of iterations for the Gibbs
//' @param adjacency A binary matrix indicating connected units
//' @param start A vector of the initializing values of
//' @param weights A numeric vector indicating the number of neighbors for each node
//' @return A numeric vector for auxiliary covariate outcomes
//' as an element of {0,1}
//'
//'
//' @export
// [[Rcpp::export]]
IntegerVector auxVarOutcomeCpp (NumericVector beta, IntegerVector trt, arma::mat cov,
                                int N, int R, List adjacency, IntegerVector start, IntegerVector weights){


  arma::mat cov_mat = cov;
  int ncol = cov_mat.n_cols;

  IntegerVector vec = start;

  // Number of iterations
  for (int r = 0; r < R; ++r){
    // Number of people

    for (int i = 0; i < N; ++i){
      float weights_i = weights[i];
      IntegerVector whichN = adjacency[i];

      // Intercept term
      float b0 = beta[0];

      // Individual treatment term
      float b1 = beta[1]*trt[i];

      // Individual covariate terms
      float b2 = 0;
      for (int q = 2; q < 2 + ncol; ++q){
        b2 = b2 + beta[q]*cov_mat(i,q-2);
      }

      // Neighbors outcome
      NumericVector v_ss = as<NumericVector>(vec[whichN]);
      float nei_w = sum(v_ss/weights_i);
      float b3 = beta[2+ncol]*nei_w;

      // Neighbors treatment
      NumericVector t_ss = as<NumericVector>(trt[whichN]);
      float nei_w2 = sum(t_ss/weights_i);
      float b4 = beta[3+ncol]*nei_w2;

      // Neighbors covariates
      float b5 = 0;
      arma::uvec whichN_arma = adjacency[i];
      for (int q = 4 + ncol; q < 4 + ncol + ncol; ++q){
        int cov_idx = q-4-ncol;
        arma::mat cov_vec = cov_mat.cols(cov_idx,cov_idx);
        float nei_w_cov = sum(cov_vec.elem(whichN_arma)/weights_i);
        b5 = b5 + beta[q]*nei_w_cov;
      }

      float prob_Lj = R::plogis(b0 + b1 + b2 + b3 + b4 + b5 , 0, 1, 1, 0);
      vec(i) = rcpp_rbinom_one(prob_Lj);
    }
  }
  return(vec);

}


//' Run Gibbs sampler for network causal effects
//'
//' Given the specific inputs, determine auxiliary covariate
//' values using a Gibbs sampling procedure.
//'
//' @param tau A numeric vector for the intercept terms in the covariate model
//' @param rho A numeric vector for the correlation terms in the covariate model
//' @param nu A numberic vector for the neighbor terms in the covariate model
//' @param ncov An integer for the number of covariates
//' @param R An integer indicating the number of iterations for the Gibbs
//' @param N An integer indicating the size of the interconnected network
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
arma::field<arma::mat> networkGibbsOutCovCpp (NumericVector tau, NumericVector rho, NumericVector nu,
                                               int ncov, int R, int N, NumericMatrix rho_mat,
                                               List adjacency,  IntegerVector weights, arma::mat cov_mat,
                                               IntegerVector group_lengths, IntegerVector group_functions){

  int J = ncov;
  int number_of_groups = group_lengths.size();
  // Create a field class (aka list) with a pre-set amount of elements
  arma::field<arma::mat> listOfMatricesOut(R);

  // Number of iterations
  for (int r = 0; r < R; ++r){

    // Number of people
    for (int i = 0; i < N; ++i){
      float weights_i = weights[i];

      // Index covariate
      int j = 0;

      // Number of groups (of covariates)
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
          NumericVector prob_vec(group_length + 1, 1.0);
          for (int m = 0; m < group_length; ++m){
            int j_prime = j + m; //already set j to 0 above whereas it is 1 in R
            arma::vec rowVec = rho_mat(_,j_prime);
            arma::mat covVec = vectorise(cov_mat.rows(i,i));
            arma::uvec whichN = adjacency[i];
            arma::mat covAdjVec = cov_mat.cols(j_prime,j_prime);
            float weights_i = weights[i];
            float nei_w = sum(covAdjVec.elem(whichN)/weights_i);
            float prob_Lj = exp(arma::as_scalar(tau[j_prime] + dot(rowVec,covVec) + nu[j_prime]*nei_w));
            prob_vec(m) = prob_Lj;
          }

          // make the rmultinom call
          NumericVector prob_vec2 = prob_vec / sum(prob_vec);
          IntegerVector multi_out = callRMultinom(prob_vec2);
          //Rcpp::Rcout << prob_vec2 << "\n";
          // update elements in matrix via another loop
          for (int m = 0; m < group_length; ++m){
            cov_mat(i,j + m) = multi_out(m);
          }
        } else if(group_function == 1){

          // Logistic / binary case

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
    listOfMatricesOut(r) = cov_mat;
  }
  return(listOfMatricesOut);

}

//' Run Gibbs sampler one
//' @export
// [[Rcpp::export]]
NumericVector networkGibbsOut1Cpp (arma::field<arma::mat> cov_list, NumericVector beta, float p,
                                              int R, int N,
                                               List adjacency,  IntegerVector weights){

}
