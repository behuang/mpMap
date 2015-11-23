#ifndef MARKER_PATTERNS_TO_UNIQUE_VALUES_HEADER_GUARD
#define MARKER_PATTERNS_TO_UNIQUE_VALUES_HEADER_GUARD
#include "Unique.hpp"
#include <map>
#include <vector>
#include "Rcpp.h"
void markerPatternsToUniqueValues(std::map<markerEncoding, markerPatternID>& markerPatterns, std::vector<markerPatternID>& markerPatternIDs, std::vector<markerEncoding>&markerEncodings, int nFounders, int nMarkers, Rcpp::IntegerMatrix& recodedFounders, long start, long end);
SEXP markerPatternsToUniqueValuesDesign(SEXP mpcross_sexp);
#endif
