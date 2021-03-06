// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "sappp_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// logit
arma::vec logit(arma::vec p);
RcppExport SEXP _sappp_logit(SEXP pSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type p(pSEXP);
    rcpp_result_gen = Rcpp::wrap(logit(p));
    return rcpp_result_gen;
END_RCPP
}
// inv_logit
arma::vec inv_logit(arma::vec a);
RcppExport SEXP _sappp_inv_logit(SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(inv_logit(a));
    return rcpp_result_gen;
END_RCPP
}
// leslie_matrix
arma::mat leslie_matrix(arma::uvec instar_days, double surv_juv, arma::vec surv_adult, arma::vec repro);
RcppExport SEXP _sappp_leslie_matrix(SEXP instar_daysSEXP, SEXP surv_juvSEXP, SEXP surv_adultSEXP, SEXP reproSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::uvec >::type instar_days(instar_daysSEXP);
    Rcpp::traits::input_parameter< double >::type surv_juv(surv_juvSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type surv_adult(surv_adultSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type repro(reproSEXP);
    rcpp_result_gen = Rcpp::wrap(leslie_matrix(instar_days, surv_juv, surv_adult, repro));
    return rcpp_result_gen;
END_RCPP
}

RcppExport SEXP _rcpp_module_boot_sappp_module();

static const R_CallMethodDef CallEntries[] = {
    {"_sappp_logit", (DL_FUNC) &_sappp_logit, 1},
    {"_sappp_inv_logit", (DL_FUNC) &_sappp_inv_logit, 1},
    {"_sappp_leslie_matrix", (DL_FUNC) &_sappp_leslie_matrix, 4},
    {"_rcpp_module_boot_sappp_module", (DL_FUNC) &_rcpp_module_boot_sappp_module, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_sappp(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
