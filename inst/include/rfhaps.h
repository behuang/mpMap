#ifndef _RFHAPS_H
#define _RFHAPS_H

#include <Rcpp.h>
#include "Unique.hpp"
#include "rfhaps_common.h"

//Return the row number within the pedigree which contains the given ID. We do this first by checking the row number the same as the ID, and then if that's not correct by just checking the whole pedigree.
inline int findIDInPedigree(int id, Rcpp::IntegerMatrix& pedigree)
{
	if(pedigree(id, 0) == id)
	{
		return id;
	}
	else
	{
		int* i = std::find(&(pedigree(0, 0)), &(pedigree(0, 0)) + pedigree.nrow(), id);
		if(i == &(pedigree(0, 0)) + pedigree.nrow()) return -1;
		return i - &(pedigree(0, 0));
	}
}
RcppExport SEXP rfhaps(SEXP RmpcrossList, SEXP recombinationFractions, SEXP marker1Range, SEXP marker2Range, SEXP RlineWeights, SEXP RuseGPU, SEXP RdeviceNum);


#endif

