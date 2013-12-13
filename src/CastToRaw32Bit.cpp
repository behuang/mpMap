#include "CastToRaw32Bit.h"

SEXP CastToRaw32Bit(SEXP input)
{
	if(TYPEOF(input) != INTSXP)
	{
		Rf_error("Input to CastToRaw32Bit must have type raw");
		return R_NilValue;
	}
	Rcpp::IntegerVector xa(input);
	Rcpp::RawVector ret(xa.size()*4);
	Rbyte* dest = &(ret[0]);
	for(Rcpp::IntegerVector::iterator i = xa.begin(); i != xa.end(); i++)
	{
		int ii = *i;
		*dest = ((ii &0xff000000) >> 24);
		dest++;
		*dest = ((ii & 0xff0000) >> 16);
		dest++;
		*dest = ((ii & 0xff00) >> 8);
		dest++;
		*dest = ((ii & 0xff) >> 0);
		dest++;
	}
	return(ret);
}
