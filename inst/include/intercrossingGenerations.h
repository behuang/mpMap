#ifndef INTERCROSSING_GENERATIONS_HEADER
#define INTERCROSSING_GENERATIONS_HEADER
#include <Rcpp.h>
//Work out the number of intercrossing generations for each individual in the final population
bool intercrossingGenerations(Rcpp::DataFrame& pedigree, int nFounders, Rcpp::IntegerVector& mpcrossID, std::vector<int>& output);
#endif
