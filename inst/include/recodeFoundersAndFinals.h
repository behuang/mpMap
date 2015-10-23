#ifndef RECODE_FOUNDERS_AND_FINALS_HEADER_GUARD
#define RECODE_FOUNDERS_AND_FINALS_HEADER_GUARD
#include "Rcpp.h"
void recodeFoundersAndFinals(Rcpp::IntegerMatrix& recodedFounders, Rcpp::IntegerMatrix& recodedFinals, Rcpp::IntegerMatrix& foundersMatrix, Rcpp::IntegerMatrix& finalsMatrix, unsigned int& maxAlleles);
#endif
