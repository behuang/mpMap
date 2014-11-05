#ifndef _RFHAPS_SPECIFIC_DESIGN_H
#define _RFHAPS_SPECIFIC_DESIGN_H
#include <Rcpp.h>
#include "Unique.hpp"
#include "sharedArray.hpp"

bool rfhapsSpecificDesign(SEXP finals, SEXP founders, SEXP pedigree_, SEXP id, SEXP fid_, SEXP recombinationFractions, long marker1Start, long marker1End, long marker2Start, long marker2End, std::vector<double>& lineWeights, SEXP RuseGPU, SEXP RdeviceNum, std::shared_ptr<double> result, std::string* error);
#endif

