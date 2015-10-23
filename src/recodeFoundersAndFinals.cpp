#include "recodeFoundersAndFinals.h"
/*
	Re-code the founder and final marker genotypes so that they always start at 0 and go up to n-1 where n is the number of distinct marker alleles at that particular marker. The maximum number of alleles across all the markers is recorded and output. 
*/
void recodeFoundersAndFinals(Rcpp::IntegerMatrix& recodedFounders, Rcpp::IntegerMatrix& recodedFinals, Rcpp::IntegerMatrix& foundersMatrix, Rcpp::IntegerMatrix& finalsMatrix, unsigned int& maxAlleles)
{
	long nMarkers = finalsMatrix.ncol();
	long nFounders = foundersMatrix.nrow(), nFinals = finalsMatrix.nrow();
	for(long markerCounter = 0; markerCounter < nMarkers; markerCounter++)
	{
		//map to hold the translation from old values to recoded values (0 - (n-1))
		std::map<int, int> translations;
		std::map<int, int>::iterator lookup;
		for(long founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			int oldValue = foundersMatrix(founderCounter, markerCounter);
			lookup = translations.find(oldValue);
			if(lookup == translations.end()) 
			{
				recodedFounders(founderCounter, markerCounter) = translations.size();
				translations.insert(std::make_pair(oldValue, translations.size()));
			}
			else
			{
				recodedFounders(founderCounter, markerCounter) = lookup->second;
			}
		}
		if(translations.size() > maxAlleles) maxAlleles = translations.size();
		//now translate founders
		for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
		{
			if(finalsMatrix(finalCounter, markerCounter) != NA_INTEGER)
			{
				recodedFinals(finalCounter, markerCounter) = translations.find(finalsMatrix(finalCounter, markerCounter))->second;
			}
			else recodedFinals(finalCounter, markerCounter) = NA_INTEGER;
		}
	}
}

