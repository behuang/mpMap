#ifndef _ALLOWABLE_MARKER_PATTERNS_HEADER_GUARD
#define _ALLOWABLE_MARKER_PATTERNS_HEADER_GUARD
#include "rfhaps.h"
void getAllowableMarkerPatterns(bool* allowableMarkerPatternsPtr, bool* allowableMarkerPatternsIRIPPtr, std::map<markerEncoding, markerPatternID>& markerPatterns, unsigned int nFounders);
#endif
