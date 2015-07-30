#ifndef VALIDATE_MPCROSS_HEADER
#define VALIDATE_MPCROSS_HEADER
#include <Rcpp.h>
#include <string>
bool validatePedigree(Rcpp::RObject pedigree_, Rcpp::DataFrame& pedigree, std::string& error);
bool validateMPCross(Rcpp::RObject mpcross_, int& nFounders, std::string& error, bool checkPedigree = true, bool checkRF = false, bool checkLG = false, bool checkFID = true);
extern "C"
{
	bool validateMPCrossNoError(SEXP mpcross_, int& nFounders, bool checkPedigree = true, bool checkRF = false, bool checkLG = false, bool checkFID = true);
}
RcppExport SEXP validate(SEXP object, SEXP checkPedigree, SEXP checkRF, SEXP checkLG, SEXP checkFID);
#endif