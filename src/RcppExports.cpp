// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// oneMultinomCall
IntegerVector oneMultinomCall(NumericVector probs);
RcppExport SEXP _autognet_oneMultinomCall(SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(oneMultinomCall(probs));
    return rcpp_result_gen;
END_RCPP
}
// oneMultinomC
IntegerVector oneMultinomC(NumericVector probs);
RcppExport SEXP _autognet_oneMultinomC(SEXP probsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type probs(probsSEXP);
    rcpp_result_gen = Rcpp::wrap(oneMultinomC(probs));
    return rcpp_result_gen;
END_RCPP
}
// callRMultinom
IntegerVector callRMultinom(NumericVector x);
RcppExport SEXP _autognet_callRMultinom(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(callRMultinom(x));
    return rcpp_result_gen;
END_RCPP
}
// auxVarCpp
arma::mat auxVarCpp(NumericVector tau, NumericVector rho, NumericVector nu, int N, int R, int J, NumericMatrix rho_mat, List adjacency, arma::mat cov_i, IntegerVector weights, IntegerVector group_lengths, IntegerVector group_functions);
RcppExport SEXP _autognet_auxVarCpp(SEXP tauSEXP, SEXP rhoSEXP, SEXP nuSEXP, SEXP NSEXP, SEXP RSEXP, SEXP JSEXP, SEXP rho_matSEXP, SEXP adjacencySEXP, SEXP cov_iSEXP, SEXP weightsSEXP, SEXP group_lengthsSEXP, SEXP group_functionsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type J(JSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rho_mat(rho_matSEXP);
    Rcpp::traits::input_parameter< List >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_i(cov_iSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group_lengths(group_lengthsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group_functions(group_functionsSEXP);
    rcpp_result_gen = Rcpp::wrap(auxVarCpp(tau, rho, nu, N, R, J, rho_mat, adjacency, cov_i, weights, group_lengths, group_functions));
    return rcpp_result_gen;
END_RCPP
}
// auxVarOutcomeCpp
IntegerVector auxVarOutcomeCpp(NumericVector beta, IntegerVector trt, arma::mat cov, int N, int R, List adjacency, IntegerVector start, IntegerVector weights);
RcppExport SEXP _autognet_auxVarOutcomeCpp(SEXP betaSEXP, SEXP trtSEXP, SEXP covSEXP, SEXP NSEXP, SEXP RSEXP, SEXP adjacencySEXP, SEXP startSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type trt(trtSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov(covSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< List >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type start(startSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(auxVarOutcomeCpp(beta, trt, cov, N, R, adjacency, start, weights));
    return rcpp_result_gen;
END_RCPP
}
// networkGibbsOutCovCpp
List networkGibbsOutCovCpp(NumericVector tau, NumericVector rho, NumericMatrix nu, int ncov, int R, int N, int burnin, NumericMatrix rho_mat, List adjacency, IntegerVector weights, arma::mat cov_mat, IntegerVector group_lengths, IntegerVector group_functions, int additional_nu);
RcppExport SEXP _autognet_networkGibbsOutCovCpp(SEXP tauSEXP, SEXP rhoSEXP, SEXP nuSEXP, SEXP ncovSEXP, SEXP RSEXP, SEXP NSEXP, SEXP burninSEXP, SEXP rho_matSEXP, SEXP adjacencySEXP, SEXP weightsSEXP, SEXP cov_matSEXP, SEXP group_lengthsSEXP, SEXP group_functionsSEXP, SEXP additional_nuSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type nu(nuSEXP);
    Rcpp::traits::input_parameter< int >::type ncov(ncovSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type rho_mat(rho_matSEXP);
    Rcpp::traits::input_parameter< List >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cov_mat(cov_matSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group_lengths(group_lengthsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type group_functions(group_functionsSEXP);
    Rcpp::traits::input_parameter< int >::type additional_nu(additional_nuSEXP);
    rcpp_result_gen = Rcpp::wrap(networkGibbsOutCovCpp(tau, rho, nu, ncov, R, N, burnin, rho_mat, adjacency, weights, cov_mat, group_lengths, group_functions, additional_nu));
    return rcpp_result_gen;
END_RCPP
}
// networkGibbsOuts1Cpp
NumericVector networkGibbsOuts1Cpp(List cov_list, NumericVector beta, float p, IntegerVector a_fixed, NumericVector dynamic_coef_vec, int dynamic_single_edge, int ncov, int R, int N, List adjacency, IntegerVector weights, IntegerVector treated_indicator, int burnin, int average, NumericVector p_vec);
RcppExport SEXP _autognet_networkGibbsOuts1Cpp(SEXP cov_listSEXP, SEXP betaSEXP, SEXP pSEXP, SEXP a_fixedSEXP, SEXP dynamic_coef_vecSEXP, SEXP dynamic_single_edgeSEXP, SEXP ncovSEXP, SEXP RSEXP, SEXP NSEXP, SEXP adjacencySEXP, SEXP weightsSEXP, SEXP treated_indicatorSEXP, SEXP burninSEXP, SEXP averageSEXP, SEXP p_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type cov_list(cov_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< float >::type p(pSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type a_fixed(a_fixedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dynamic_coef_vec(dynamic_coef_vecSEXP);
    Rcpp::traits::input_parameter< int >::type dynamic_single_edge(dynamic_single_edgeSEXP);
    Rcpp::traits::input_parameter< int >::type ncov(ncovSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< List >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type treated_indicator(treated_indicatorSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type average(averageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p_vec(p_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(networkGibbsOuts1Cpp(cov_list, beta, p, a_fixed, dynamic_coef_vec, dynamic_single_edge, ncov, R, N, adjacency, weights, treated_indicator, burnin, average, p_vec));
    return rcpp_result_gen;
END_RCPP
}
// networkGibbsOuts2Cpp
NumericVector networkGibbsOuts2Cpp(List cov_list, NumericVector beta, float p, IntegerVector a_fixed, NumericVector dynamic_coef_vec, int dynamic_single_edge, int ncov, int R, int N, List adjacency, IntegerVector weights, IntegerVector treated_indicator, IntegerVector subset, float treatment_value, int burnin, int average, NumericVector p_vec);
RcppExport SEXP _autognet_networkGibbsOuts2Cpp(SEXP cov_listSEXP, SEXP betaSEXP, SEXP pSEXP, SEXP a_fixedSEXP, SEXP dynamic_coef_vecSEXP, SEXP dynamic_single_edgeSEXP, SEXP ncovSEXP, SEXP RSEXP, SEXP NSEXP, SEXP adjacencySEXP, SEXP weightsSEXP, SEXP treated_indicatorSEXP, SEXP subsetSEXP, SEXP treatment_valueSEXP, SEXP burninSEXP, SEXP averageSEXP, SEXP p_vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< List >::type cov_list(cov_listSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< float >::type p(pSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type a_fixed(a_fixedSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type dynamic_coef_vec(dynamic_coef_vecSEXP);
    Rcpp::traits::input_parameter< int >::type dynamic_single_edge(dynamic_single_edgeSEXP);
    Rcpp::traits::input_parameter< int >::type ncov(ncovSEXP);
    Rcpp::traits::input_parameter< int >::type R(RSEXP);
    Rcpp::traits::input_parameter< int >::type N(NSEXP);
    Rcpp::traits::input_parameter< List >::type adjacency(adjacencySEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type treated_indicator(treated_indicatorSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type subset(subsetSEXP);
    Rcpp::traits::input_parameter< float >::type treatment_value(treatment_valueSEXP);
    Rcpp::traits::input_parameter< int >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< int >::type average(averageSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type p_vec(p_vecSEXP);
    rcpp_result_gen = Rcpp::wrap(networkGibbsOuts2Cpp(cov_list, beta, p, a_fixed, dynamic_coef_vec, dynamic_single_edge, ncov, R, N, adjacency, weights, treated_indicator, subset, treatment_value, burnin, average, p_vec));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_autognet_oneMultinomCall", (DL_FUNC) &_autognet_oneMultinomCall, 1},
    {"_autognet_oneMultinomC", (DL_FUNC) &_autognet_oneMultinomC, 1},
    {"_autognet_callRMultinom", (DL_FUNC) &_autognet_callRMultinom, 1},
    {"_autognet_auxVarCpp", (DL_FUNC) &_autognet_auxVarCpp, 12},
    {"_autognet_auxVarOutcomeCpp", (DL_FUNC) &_autognet_auxVarOutcomeCpp, 8},
    {"_autognet_networkGibbsOutCovCpp", (DL_FUNC) &_autognet_networkGibbsOutCovCpp, 14},
    {"_autognet_networkGibbsOuts1Cpp", (DL_FUNC) &_autognet_networkGibbsOuts1Cpp, 15},
    {"_autognet_networkGibbsOuts2Cpp", (DL_FUNC) &_autognet_networkGibbsOuts2Cpp, 17},
    {NULL, NULL, 0}
};

RcppExport void R_init_autognet(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
