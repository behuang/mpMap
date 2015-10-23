#include "markerPatternsToUniqueValues.h"
#include "recodeFoundersAndFinals.h"
#include "validateMPCross.h"
bool sortMarkerPatternIDs(const std::pair<long, markerPatternID>& first, const std::pair<long, markerPatternID>& second)
{
	return (int)first.second < (int)second.second;
}
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
	std::vector<markerPatternID> markerPatternIDs;
	//vector with entry i containing an encoding of the marker segregation pattern for a marker with markerPatternID i
	std::vector<markerEncoding> markerEncodings;
	markerPatternsToUniqueValues(markerPatterns, markerPatternIDs, markerEncodings, nFounders, nMarkers, recodedFounders);

	//We're going to make a vector of pairs, containing markerPatternIDs and its index
	std::vector<std::pair<long, markerPatternID> > sortedMarkerPatternIDs;
	sortedMarkerPatternIDs.reserve(markerPatternIDs.size());
	
	for(std::vector<markerPatternID>::iterator i = markerPatternIDs.begin(); i != markerPatternIDs.end(); i++)
	{
		std::size_t index = std::distance(markerPatternIDs.begin(), i);
		sortedMarkerPatternIDs.push_back(std::make_pair((long)index, *i));
	}
	std::sort(sortedMarkerPatternIDs.begin(), sortedMarkerPatternIDs.end(), sortMarkerPatternIDs);

	//Now get out just the indices
	Rcpp::IntegerVector permutation(nMarkers);
	for(std::vector<std::pair<long, markerPatternID> >::iterator i = sortedMarkerPatternIDs.begin(); i != sortedMarkerPatternIDs.end(); i++)
	{
		std::size_t index = std::distance(sortedMarkerPatternIDs.begin(), i);
		permutation[(int)index] = i->first+1;
	}
	
	return permutation;
END_RCPP
}
void markerPatternsToUniqueValues(std::map<markerEncoding, markerPatternID>& markerPatterns, std::vector<markerPatternID>& markerPatternIDs, std::vector<markerEncoding>&markerEncodings, int nFounders, int nMarkers, Rcpp::IntegerMatrix& recodedFounders)
{
	for(long markerCounter = 0; markerCounter < nMarkers; markerCounter++)
	{
		int encodedMarkerPattern = 0;
		for(long founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			encodedMarkerPattern += (recodedFounders(founderCounter, markerCounter) << 3*founderCounter);
		}
		if(markerPatterns.find(encodedMarkerPattern) == markerPatterns.end())
		{
			markerPatternIDs.push_back(markerPatterns.size());
			markerPatterns.insert(std::make_pair(encodedMarkerPattern, markerPatterns.size()));
			markerEncodings.push_back(encodedMarkerPattern);
		}
		else
		{
			markerPatternIDs.push_back(markerPatterns.find(encodedMarkerPattern)->second);
		}
	}
}

