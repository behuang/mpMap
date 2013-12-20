#ifndef GENERATE_GENOTYPES_HEADER_GUARD
#define GENERATE_GENOTYPES_HEADER_GUARD
#include <Rcpp.h>
RcppExport SEXP generateGenotypes(SEXP RrecombinationDistances, SEXP Rpedigree, SEXP RtransMarkerIndices, SEXP RtransFounder);
#endif