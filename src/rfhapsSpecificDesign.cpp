#include "rfhapsSpecificDesign.h"
#include <math.h>
#include <sstream>
#include "sharedArray.hpp"
#include "rfhaps_cpu.h"
#include "rfhaps_gpu.h"
#include <limits>
#include "getFunnelCPU.h"
#include "validateMPCross.h"
#include "intercrossingGenerations.h"
#include "orderFunnels.h"
#include "allowableMarkerPatterns.h"
#include <array>
#include "findIDInPedigree.h"
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
/*
	AIC lines don't have a unique funnel. But if we go back far enough we'll get a bunch of individuals, each of which DOES have a unique funnel. This function gets out the rows in the pedigree corresponding to ancestors of the line with the given id. 
	These ancestors have unique funnels. 
*/
void getAICParentsWithFunnels(Rcpp::DataFrame& pedigreeDataFrame, long id, int nIntercrossingGenerations, std::vector<long>& individualsToCheckFunnels)
{
	Rcpp::IntegerVector ids = pedigreeDataFrame("id");
	Rcpp::IntegerVector male = pedigreeDataFrame("Male");
	Rcpp::IntegerVector female = pedigreeDataFrame("Female");
	//The lines that we currently need to check goes in individualsToCheckFunnels
	individualsToCheckFunnels.clear();
	//convert initial ID into a row in the pedigree
	int currentPedigreeRow = findIDInPedigree(id, pedigreeDataFrame);
	//step back through all the selfing generations, if they're explicitly listed
	while(male(currentPedigreeRow) == female(currentPedigreeRow))
	{
		currentPedigreeRow = findIDInPedigree(male(currentPedigreeRow), pedigreeDataFrame);
	}
	individualsToCheckFunnels.push_back(currentPedigreeRow);

	//step back through the intercrossing generations
	std::vector<long> nextGenerationToCheck;
	for(;nIntercrossingGenerations > 0; nIntercrossingGenerations--)
	{
		nextGenerationToCheck.clear();
		//For everything that currently needs to be checked, go up one generation. 
		for(std::vector<long>::iterator i = individualsToCheckFunnels.begin(); i != individualsToCheckFunnels.end(); i++)
		{
			nextGenerationToCheck.push_back(findIDInPedigree(male(*i), pedigreeDataFrame));
			nextGenerationToCheck.push_back(findIDInPedigree(female(*i), pedigreeDataFrame));
		}
		//swap vectors
		individualsToCheckFunnels.swap(nextGenerationToCheck);
	}
	//At the end the vector contains PEDIGREE ROWS not IDs, so change to IDs
	nextGenerationToCheck.clear();
	for(std::vector<long>::iterator i = individualsToCheckFunnels.begin(); i != individualsToCheckFunnels.end(); i++)
	{
		nextGenerationToCheck.push_back(ids(*i));
	}
	individualsToCheckFunnels.swap(nextGenerationToCheck);
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
void funnelsToUniqueValues(std::map<funnelEncoding, funnelID>& funnelTranslation, std::vector<funnelID>& funnelIDs, std::vector<funnelEncoding>& funnelEncodings, std::vector<intArray8>& allFunnels, int nFounders)
{
	for(std::vector<intArray8>::iterator i = allFunnels.begin(); i != allFunnels.end(); i++)
	{
		intArray8 funnel = *i;
		int encoded = 0;
		for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
		{
			encoded += ((funnel.val[founderCounter]-1) << (3 * founderCounter));
		}
		if(funnelTranslation.find(encoded) == funnelTranslation.end())
		{
			funnelIDs.push_back(funnelTranslation.size());
			funnelTranslation.insert(std::make_pair(encoded, funnelTranslation.size()));
			funnelEncodings.push_back(encoded);
		}
		else
		{
			funnelIDs.push_back(funnelTranslation.find(encoded)->second);
		}
	}
}
bool rfhapsSpecificDesign(SEXP finals, SEXP founders, SEXP pedigree_, SEXP id, SEXP fid_, SEXP recombinationFractions, long marker1Start, long marker1End, long marker2Start, long marker2End, std::vector<double>& lineWeights, SEXP RuseGPU, SEXP RdeviceNum, std::shared_ptr<double> result, std::string* error)
{
	//Convert SEXP to R data types. No error checking here because it's been checked previously in function rfhaps
	Rcpp::IntegerMatrix foundersMatrix(founders);
	Rcpp::IntegerMatrix finalsMatrix(finals);
	Rcpp::DataFrame pedigreeDataFrame(pedigree_);
	Rcpp::NumericVector recombinationVector(recombinationFractions);
	Rcpp::IntegerVector useGPU_(RuseGPU);
	Rcpp::IntegerVector deviceNum_(RdeviceNum);
	Rcpp::IntegerVector IDVector(id);
	Rcpp::IntegerVector fidVector(fid_);

	Rcpp::List finalsDimnames = finalsMatrix.attr("dimnames");
	Rcpp::CharacterVector allMarkerNames = finalsDimnames[1];
	
	long nMarkers = finalsMatrix.ncol();

	bool useGPU = useGPU_(0);
	int deviceNum = deviceNum_(0);


	long nFounders = foundersMatrix.nrow(), nFinals = finalsMatrix.nrow();
	if(IDVector.size() != nFinals)
	{
		*error = "Entry id had the wrong number of entries";
		return false;
	}
	//work out if there are any AI generations. Computation simplifies if there aren't any
	std::vector<int> nIntercrossingGenerations;
	nIntercrossingGenerations.resize(nFinals, 0);
	//Returns false if the pedigree was invalid / had wrong form
	bool pedigreeOK = intercrossingGenerations(pedigreeDataFrame, nFounders, IDVector, nIntercrossingGenerations);
	if(!pedigreeOK)
	{
		*error = "There was a problem with the input pedigree, terminating";
		return false;
	}
	bool hasAI = false;
	for(std::vector<int>::iterator i = nIntercrossingGenerations.begin(); i != nIntercrossingGenerations.end(); i++)
	{
		if(*i > 0) 
		{
			hasAI = true; 
			break;
		}
	}
	if(nFounders != 2 && nFounders != 4 && nFounders != 8)
	{
		*error = "Number of founders must be 2, 4 or 8";
		return false;
	}
	if(fidVector.length() != nFounders)
	{
		*error = "Input mpcross$fid had the wrong length";
		return false;
	}
	std::vector<int> fid = Rcpp::as<std::vector<int> >(fidVector);
	//Check the funnels.
	//We have a vector of individuals to check, because we also want to check the parents of the AIC lines, which only have well-defined funnels further up the pedigree, so we need to check several at a time in that case
	std::vector<long> individualsToCheckFunnels;
	//The number of warning messages that have been printed
	int messages = 0;
	//All the funnels present in the non-intercrossing lines
	std::vector<intArray8> allFunnels;
	for(long individualCounter = 0; individualCounter < nFinals; individualCounter++)
	{
		individualsToCheckFunnels.clear();
		if(nIntercrossingGenerations[individualCounter] == 0)
		{
			individualsToCheckFunnels.push_back(IDVector(individualCounter));
		}
		else
		{
			getAICParentsWithFunnels(pedigreeDataFrame, IDVector(individualCounter), nIntercrossingGenerations[individualCounter], individualsToCheckFunnels);
		}
		//Now we know the lines for which we need to check the funnels from the pedigree (note: We don't necessarily have genotype data for all of these, it's purely a pedigree check)
		std::vector<int> representedFounders;
		//Fixed length arrays to store funnels. If we have less than 8 founders then part of this is garbage and we don't use that bit....
		intArray8 funnel, copiedFunnel;
		for(std::vector<long>::iterator i = individualsToCheckFunnels.begin(); i != individualsToCheckFunnels.end(); i++)
		{
			if(!getFunnel(*i, pedigreeDataFrame, fid, 0, &(funnel.val[0]), pedigreeDataFrame.nrows(), nFounders))
			{
				std::stringstream ss;
				ss << "Problem with pedigree for individual %d\n" << (individualCounter+1);
				*error = ss.str();
				return false;
			}
			//insert these founders into the vector containing all the represented founders
			representedFounders.insert(representedFounders.end(), &(funnel.val[0]), &(funnel.val[0]) + nFounders);
			//Copy the funnel 
			memcpy(&copiedFunnel, &funnel, sizeof(intArray8));
			std::sort(&(copiedFunnel.val[0]), &(copiedFunnel.val[0]) + nFounders);
			if(std::unique(&(copiedFunnel.val[0]), &(copiedFunnel.val[0]) + nFounders) != &(copiedFunnel.val[0]) + nFounders && messages < 6)
			{
				//Not having all the founders represented is only a warning, by itself....
				Rprintf("Warning: Funnel for individual %d contained founders {%d, ", *i, funnel.val[0]);
				if(nFounders == 2)
				{
					Rprintf("%d}", funnel.val[1]);
				}
				else if(nFounders == 4)
				{
					Rprintf("%d, %d, %d}", funnel.val[1], funnel.val[2], funnel.val[3]);
				}
				else
				{
					Rprintf("%d, %d, %d, %d, %d, %d, %d}", funnel.val[1], funnel.val[2], funnel.val[3], funnel.val[4], funnel.val[5], funnel.val[6], funnel.val[7]);
				}
				Rprintf(". Did you intend to use all %d founders?\n", nFounders);
				messages++;
			}
			if(messages == 6) Rprintf("Suppressing further warnings\n");
		}
		//remove duplicates in representedFounders
		std::sort(representedFounders.begin(), representedFounders.end());
		representedFounders.erase(std::unique(representedFounders.begin(), representedFounders.end()), representedFounders.end());
		//Not having all the founders in the input funnels is more serious if it causes the observed marker data to be impossible. So check for this.
		for(int markerCounter = 0; markerCounter < nMarkers; markerCounter++)
		{
			bool okMarker = false;
			//If observed value is an NA then than's ok, continue
			int value = finalsMatrix(individualCounter, markerCounter);
			if(value == NA_INTEGER) continue;

			for(std::vector<int>::iterator founderIterator = representedFounders.begin(); founderIterator != representedFounders.end(); founderIterator++)
			{
				if(finalsMatrix(individualCounter, markerCounter) == foundersMatrix((*founderIterator)-1, markerCounter))
				{
					okMarker = true;
					break;
				}
			}
			if(!okMarker)
			{
				std::stringstream ss;
				ss << "Error: Data for marker number " << markerCounter + 1 << " is impossible for individual " << individualCounter + 1 << " with given pedigree\n";
				*error = ss.str();
				return false;
			}
		}
		//In this case individualsToCheckFunnels contains one element => getFunnel was only called once => we can reuse the funnel variable
		if(nIntercrossingGenerations[individualCounter] == 0)
		{
			orderFunnel(&(funnel.val[0]), nFounders);
			allFunnels.push_back(funnel);
		}
	}
	//re-code the founder and final marker genotypes so that they always start at 0 and go up to n-1 where n is the number of distinct marker alleles
	//We do this to make it easier to identify markers with identical segregation patterns. recodedFounders = column major matrix
	Rcpp::IntegerMatrix recodedFounders(nFounders, nMarkers), recodedFinals(nFinals, nMarkers);
	recodedFinals.attr("dimnames") = finalsMatrix.attr("dimnames");
	unsigned int maxAlleles = 0;
	recodeFoundersAndFinals(recodedFounders, recodedFinals, foundersMatrix, finalsMatrix, maxAlleles);
	if(maxAlleles > 8) 
	{
		*error = "Internal error - Cannot have more than eight alleles per marker";
		return false;
	}
	//map containing encodings of all the haplotypes, and an associated unique index (we can't use the haplotype encoding as an index because they'll jump around a lot and could be quite big values sometimes I think). Unique indices are guaranteed to be contiguous numbers starting from 0 - So [0, markerPatterns.size()]. 
	std::map<markerEncoding, markerPatternID> markerPatterns;
	//A vector where entry i contains the markerPatternID identifying the segregation pattern of marker number i. Contains one entry per marker.
	std::vector<markerPatternID> markerPatternIDs;
	//vector with entry i containing an encoding of the marker segregation pattern for a marker with markerPatternID i
	std::vector<markerEncoding> markerEncodings;
	markerPatternsToUniqueValues(markerPatterns, markerPatternIDs, markerEncodings, nFounders, nMarkers, recodedFounders);
	
	//map containing encodings of the funnels involved in the experiment (as key), and an associated unique index (again, using the encoded values directly is no good because they'll be all over the place). Unique indices are contiguous again.
	std::map<funnelEncoding, funnelID> funnelTranslation;
	//vector giving the funnel ID for each individual
	std::vector<funnelID> funnelIDs;
	//vector giving the encoded value for each individual
	std::vector<funnelEncoding> funnelEncodings;
	funnelsToUniqueValues(funnelTranslation, funnelIDs, funnelEncodings, allFunnels, nFounders);
	
	//Construct boolean matrix where rows and columns represent marker segregation patterns, and the boolean values refer to whether or not that pair of marker segregation patterns can be used to estimate recombination fractions
	//A pair can be unaccetable as all parameters lead to the same probability model (complete unidentifiability) or there are pairs of parameters that lead to the same probability model (We will be able to estimate the "best pair", but get no further).
	sharedArray<bool> allowableMarkerPatterns(new bool[markerPatterns.size() * markerPatterns.size()]);
	bool* allowableMarkerPatternsPtr = allowableMarkerPatterns.get();

	getAllowableMarkerPatterns(allowableMarkerPatternsPtr, markerPatterns, nFounders);
	bool resultOK;

	rfhaps_cpu_args cpu_args(nIntercrossingGenerations, markerPatternIDs, lineWeights, markerEncodings, funnelIDs, funnelEncodings);
	cpu_args.foundersMatrix = recodedFounders;
	cpu_args.finalsMatrix = recodedFinals;
	cpu_args.pedigreeMatrix = pedigreeDataFrame;
	cpu_args.recombinationVector = recombinationVector;
	cpu_args.IDVector = IDVector;
	cpu_args.fidVector = fidVector;
	cpu_args.nFounders = nFounders;
	cpu_args.hasAI = hasAI;
	cpu_args.result = result.get();
	cpu_args.marker1Start = marker1Start;
	cpu_args.marker2Start = marker2Start;
	cpu_args.marker1End = marker1End;
	cpu_args.marker2End = marker2End;
	cpu_args.allowableMarkerPatterns = allowableMarkerPatterns;
	cpu_args.maxMarkerAlleles = maxAlleles;
#ifndef HAS_CUDA
	if(useGPU)
	{
		Rprintf("Package was not compiled with GPU support, falling back to CPU\n");
	}
	resultOK = rfhaps_cpu(cpu_args);
#else
	if(!useGPU)
	{
		resultOK = rfhaps_cpu(cpu_args);
	}
	else
	{
		std::vector<std::string> design;
		if(pedigreeDataFrame.length() == 5)
		{
			design = Rcpp::as<std::vector<std::string> >(pedigreeDataFrame("Design"));
		}
		pedigreeColumns pedigree(&(Rcpp::as<Rcpp::IntegerVector>(pedigreeDataFrame("id"))(0)), &(Rcpp::as<Rcpp::IntegerVector>(pedigreeDataFrame("Male"))(0)), &(Rcpp::as<Rcpp::IntegerVector>(pedigreeDataFrame("Female"))(0)), &(Rcpp::as<Rcpp::IntegerVector>(pedigreeDataFrame("Observed"))(0)), design);
		rfhaps_gpu_args args(pedigree, markerPatternIDs, fid, lineWeights, markerEncodings, funnelIDs, funnelEncodings, allFunnels);
		args.founders = &(recodedFounders(0, 0));
		args.finals = &(recodedFinals(0, 0));
		args.recombination = &(recombinationVector(0));
		args.nIntercrossing = &(nIntercrossingGenerations[0]);
		args.nFounders = nFounders;
		args.nFinals = nFinals;
		args.nRecomb = recombinationVector.length();
		args.hasAI = hasAI;
		args.output = result.get();
		args.marker1Start = marker1Start;
		args.marker1End = marker1End;
		args.marker2Start = marker2Start;
		args.marker2End = marker2End;
		args.nPedigreeRows = pedigreeDataFrame.nrows();
		args.allowableMarkerPatterns = allowableMarkerPatternsPtr;
		args.IDs = &(IDVector(0));
		args.deviceNum = deviceNum;
		resultOK = rfhaps_gpu(args);
	}
#endif
	if(!resultOK) return false;
	return true;
}
