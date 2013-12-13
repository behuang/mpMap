#include <Rcpp.h>
#include "rfhaps.h"
#include "CastToRaw32Bit.h"
#include "impute.h"
#include "getAllFunnels.h"
extern "C"
{
	R_CallMethodDef callMethods[] = 
	{
		{"rfhaps", (DL_FUNC)&rfhaps, 6},
		{"CastToRaw32Bit", (DL_FUNC)CastToRaw32Bit, 1},
		{"impute", (DL_FUNC)impute, 1},
		{"getAllFunnels", (DL_FUNC)getAllFunnels, 1},
		{NULL, NULL, 0}
	};
	RcppExport void R_init_EightWay9k(DllInfo *info)
	{
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	}
}
