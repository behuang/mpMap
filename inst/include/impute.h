#ifndef _EightWay9k_IMPUTE_H
#define _EightWay9k_IMPUTE_H

#include <Rcpp.h>
extern "C"
{
	bool imputeInternal(double* theta, double* lod, double* lkhd, int nMarkers, int* groups, char* error, int errorLength);
}
RcppExport SEXP impute(SEXP mpcross);
#endif

