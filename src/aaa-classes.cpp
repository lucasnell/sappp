#include <RcppArmadillo.h> // arma namespace
#include <pcg/pcg_random.hpp> // pcg prng
#include <vector>          // vector class
#include <cmath>           // log, exp
#include <random>          // normal distribution
#include <cstdint>         // integer types

#include "sappp_types.h"   // integer types
#include "cpp_classes.h"   // all these classes

using namespace Rcpp;


RCPP_EXPOSED_CLASS(SimSummary)
RCPP_EXPOSED_CLASS(SimPatches)

    

RCPP_MODULE(sappp_module) {
    
    class_<SimSummary>("SimSummary")
        .field("aphids", &SimSummary::aphids, "density of unparasitized aphids")
        .field("parasit", &SimSummary::parasit, 
            "density of parasitized, but alive, aphids")
        .field("mummies", &SimSummary::mummies, "density of mummies")
        .field("wasps", &SimSummary::wasps, "density of adult wasps")
        .method("show", &SimSummary::show, "summarize SimSummary object")
        .method("flatten", &SimSummary::flatten, "flatten to a single matrix")
    ;
    
    class_<SimPatches>("SimPatches")
        .constructor< std::vector<std::vector<List>>,std::vector<uint32>,std::vector<uint32> >()
        .field_readonly("aphid_names", &SimPatches::aphid_names, 
            "vector of aphid names (same for all patches)")
        .field("t", &SimPatches::t, "current time")
        .method("simulate", &SimPatches::simulate,
            "reset object and simulate a set number of time periods")
    ;
}

