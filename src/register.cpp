#include <Rcpp.h>
#include <R_ext/Rdynload.h>
#include "rfhaps.h"
#include "impute.h"
#include "getAllFunnels.h"
#include "validateMPCross.h"
#include "markerPatternsToUniqueValues.h"
extern "C"
{
	R_CallMethodDef callMethods[] = 
	{
		{"rfhaps", (DL_FUNC)&rfhaps, 7},
		{"impute", (DL_FUNC)impute, 1},
		{"getAllFunnels", (DL_FUNC)getAllFunnels, 1},
		{"validateMPCross", (DL_FUNC)validate, 5},
		{"markerPatternsToUniqueValuesDesign", (DL_FUNC)markerPatternsToUniqueValuesDesign, 1},
		{NULL, NULL, 0}
	};
	RcppExport void R_init_mpMap(DllInfo *info)
	{
		R_registerRoutines(info, NULL, callMethods, NULL, NULL);
		R_RegisterCCallable("mpMap", "imputeInternal", (DL_FUNC)imputeInternal);
		R_RegisterCCallable("mpMap", "validateMPCrossNoError", (DL_FUNC)validateMPCrossNoError);
	}
}
