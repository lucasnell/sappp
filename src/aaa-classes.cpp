#include <RcppArmadillo.h> // arma namespace
#include <sitmo.h>         // parallel rng
#include <vector>          // vector class
#include <cmath>           // log, exp
#include <random>          // normal distribution
#include <cstdint>         // integer types

#include "sap_types.h"

using namespace Rcpp;
using namespace std;


RCPP_EXPOSED_CLASS(aphid_wasp)
RCPP_MODULE(aphid_wasp_mod) {
    
    class_<aphid_wasp>("aphid_wasp")
        .constructor<List>()
        .method("show", &aphid_wasp::show)
        .field_readonly("aphid_name", &aphid_wasp::aphid_name, 
            "unique identifying name for this aphid line")
    ;
}

