#ifndef _RFHAPS_CPU_H
#define _RFHAPS_CPU_H
#include <Rcpp.h>
#include <memory>
#include "Unique.hpp"
struct rfhaps_cpu_args
{
	rfhaps_cpu_args(std::vector<int>& nIntercrossingGenerations, std::vector<markerPatternID>& markerPatternIDs, std::vector<double>& lineWeights, std::vector<markerEncoding>& markerEncodings, std::vector<funnelID>& funnelIDs, std::vector<funnelEncoding>& funnelEncodings)
	: nIntercrossingGenerations(nIntercrossingGenerations), markerPatternIDs(markerPatternIDs), lineWeights(lineWeights), markerEncodings(markerEncodings), funnelIDs(funnelIDs), funnelEncodings(funnelEncodings)
	{}
	Rcpp::IntegerMatrix foundersMatrix;
	Rcpp::IntegerMatrix finalsMatrix;
	Rcpp::DataFrame pedigreeMatrix;
	Rcpp::NumericVector recombinationVector;
	Rcpp::IntegerVector IDVector;
	Rcpp::IntegerVector fidVector;
	
	std::vector<int>& nIntercrossingGenerations;
	//A vector where entry i contains the markerPatternID identifying the segregation pattern of marker number i. Contains one entry per marker. 
	std::vector<markerPatternID>& markerPatternIDs;
	std::vector<double>& lineWeights;
	//vector with entry i containing an encoding of the marker segregation pattern for a marker with markerPatternID i
	std::vector<markerEncoding>& markerEncodings;
	//A symmetric boolean matrix, with dimensions translations.size(). This says whether segregation patterns NUMBER i and j (numbered among all segregation patterns by the translations map) can be used to estimate recombination fractions. 
	std::shared_ptr<bool> allowableMarkerPatterns;
	int nFounders;
	bool hasAI;
	int marker1Start, marker1End;
	int marker2Start, marker2End;
	//maximum number of marker alleles present
	int maxMarkerAlleles;
	//zeroed by the calling code. Must be added to, not overwritten. 
	double* result;
	std::vector<funnelID>& funnelIDs;
	std::vector<funnelEncoding>& funnelEncodings;
};
bool rfhaps_cpu(rfhaps_cpu_args& args);
#endif
