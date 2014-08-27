#ifndef VALIDATE_MPCROSS_HEADER
#define VALIDATE_MPCROSS_HEADER
#include <Rcpp.h>
#include <string>
bool validateMPCross(Rcpp::RObject mpcross_, int& nFounders, std::string& error, bool checkPedigree = true, bool checkRF = false, bool checkLG = false, bool checkFID = true);
RcppExport SEXP validate(SEXP object, SEXP checkPedigree, SEXP checkRF, SEXP checkLG, SEXP checkFID);
#endif