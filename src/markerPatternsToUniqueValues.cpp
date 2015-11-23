#include "markerPatternsToUniqueValues.h"
#include "recodeFoundersAndFinals.h"
#include "validateMPCross.h"
SEXP markerPatternsToUniqueValuesDesign(SEXP mpcross_sexp)
{
BEGIN_RCPP
	int nFounders;
	std::string error;
	bool valid = validateMPCross(mpcross_sexp, nFounders, error, true, false, false);
	if(!valid)
	{
		throw std::runtime_error(error);
	}
	Rcpp::List mpcross(mpcross_sexp);
	Rcpp::IntegerMatrix foundersMatrix = mpcross["founders"];
	Rcpp::IntegerMatrix finalsMatrix = mpcross["finals"];
	long nMarkers = finalsMatrix.ncol();
	long nFinals = finalsMatrix.nrow();	

	Rcpp::IntegerMatrix recodedFounders(nFounders, nMarkers), recodedFinals(nFinals, nMarkers);
	unsigned int maxAlleles = 0;
	recodeFoundersAndFinals(recodedFounders, recodedFinals, foundersMatrix, finalsMatrix, maxAlleles);

	std::map<markerEncoding, markerPatternID> markerPatterns;
	//A vector where entry i contains the markerPatternID identifying the segregation pattern of marker number i. Contains one entry per marker.
	std::vector<markerPatternID> markerPatternIDs(nMarkers);
	//vector with entry i containing an encoding of the marker segregation pattern for a marker with markerPatternID i
	std::vector<markerEncoding> markerEncodings;
	markerPatternsToUniqueValues(markerPatterns, markerPatternIDs, markerEncodings, nFounders, nMarkers, recodedFounders, 0, nMarkers);

	Rcpp::IntegerVector convertedMarkerPatternIDs(nMarkers);
	std::copy(markerPatternIDs.begin(), markerPatternIDs.end(), convertedMarkerPatternIDs.begin());
	return convertedMarkerPatternIDs;
END_RCPP
}
void markerPatternsToUniqueValues(std::map<markerEncoding, markerPatternID>& markerPatterns, std::vector<markerPatternID>& markerPatternIDs, std::vector<markerEncoding>&markerEncodings, int nFounders, int nMarkers, Rcpp::IntegerMatrix& recodedFounders, long start, long end)
{
	for(long markerCounter = start; markerCounter < end; markerCounter++)
	{
		int encodedMarkerPattern = 0;
		for(long founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			encodedMarkerPattern += (recodedFounders(founderCounter, markerCounter) << 3*founderCounter);
		}
		if(markerPatterns.find(encodedMarkerPattern) == markerPatterns.end())
		{
			markerPatternIDs[markerCounter] = markerPatterns.size();
			markerPatterns.insert(std::make_pair(encodedMarkerPattern, markerPatterns.size()));
			markerEncodings.push_back(encodedMarkerPattern);
		}
		else
		{
			markerPatternIDs[markerCounter] = markerPatterns.find(encodedMarkerPattern)->second;
		}
	}
}

