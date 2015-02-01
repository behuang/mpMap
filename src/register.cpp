#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include "rfhaps.h"
#include "CastToRaw32Bit.h"
#include "impute.h"
#include "getAllFunnels.h"
#include "validateMPCross.h"
extern "C"
{
	R_CallMethodDef callMethods[] = 
	{
		{"rfhaps", (DL_FUNC)&rfhaps, 7},
		{"CastToRaw32Bit", (DL_FUNC)CastToRaw32Bit, 1},
		{"impute", (DL_FUNC)impute, 1},
		{"getAllFunnels", (DL_FUNC)getAllFunnels, 1},
		{"validateMPCross", (DL_FUNC)validate, 5},
		{NULL, NULL, 0}
	};
	RcppExport void R_init_mpMap(DllInfo *info)
	{
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
		R_RegisterCCallable("mpMap", "imputeInternal", (DL_FUNC)imputeInternal);
	}
}
