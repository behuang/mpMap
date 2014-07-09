#include "rfhaps_cpu.h"
#include "rfhaps.h"
#include "sharedArray.hpp"
#include <math.h>
#include "getFunnelCPU.h"
#include <array>
//Note that if we change the mask[8][8] values of 2 to 1 we get mask4 in the first 4x4 block. 
//const int mask4[4][4] = {{0, 1, 1, 1}, {1, 0, 1, 1}, {1, 1, 0, 1}, {1, 1, 1, 0}};
const int mask[8][8] =
	{
			{0, 1, 2, 2, 2, 2, 2, 2},
			{1, 0, 2, 2, 2, 2, 2, 2},
			{2, 2, 0, 1, 2, 2, 2, 2},
			{2, 2, 1, 0, 2, 2, 2, 2},
			{2, 2, 2, 2, 0, 1, 2, 2},
			{2, 2, 2, 2, 1, 0, 2, 2},
			{2, 2, 2, 2, 2, 2, 0, 1},
			{2, 2, 2, 2, 2, 2, 1, 0}
	};
//Templated function to work out the two-point probabilities with the given recombination fraction (and number of AI generations). Templating allows the number of founders to be a compile-time constant
template<int nFounders> void genotypeProbabilitiesNoIntercross(double (&prob)[3], double recombinationFraction);
template<int nFounders> void genotypeProbabilitiesWithIntercross(double (&prob)[3], int nAIGenarations, double recombinationFraction);
//There are only really two probability values for the 4-way design, but if we put in three we can use the same mask as for the 8-way case
template<> void genotypeProbabilitiesNoIntercross<4>(double (&prob)[3], double r)
{
	prob[0] = (1-r)/(4+8*r);
	prob[1] = prob[2] = r/(4+8*r);
}
template<> void genotypeProbabilitiesNoIntercross<8>(double (&prob)[3], double r)
{
	prob[0] = (1-r)*(1-r)/(8+16*r);
	prob[1] = r*(1-r)/(8+16*r);
	prob[2] = r/(16+32*r);
}
template<> void genotypeProbabilitiesNoIntercross<2>(double (&prob)[3], double r)
{
	prob[0] = 1/(2*(1 + 2*r));
	prob[1] = r/(1 + 2 * r);
}
template<> void genotypeProbabilitiesWithIntercross<2>(double (&prob)[3], int nAIGenerations, double r)
{
	double tmp = pow(1-r, nAIGenerations - 1);
	//calculated by taking the 4-way case and setting both pairs of founders to be identical
	prob[0] = 1 - (4.0/3.0) * (1/(1 + 2 * r)) * ((1-r)*(1-r)*tmp/4 + (2*r + 1 - tmp) /16);
	prob[1] = (1 - prob[0]*2)/2;
}
template<> void genotypeProbabilitiesWithIntercross<4>(double (&prob)[3], int nAIGenerations, double r)
{
	double tmp = pow(1-r, nAIGenerations-1);
	//prob[0] = (pow(1-r, 1+nAIGenerations)/4+(2*r+1-pow(1-r, nAIGenerations-1))/16)/(1+2*r); 
	prob[0] = (tmp *(1-r)*(1-r)/4 + (2*r + 1 - tmp)/16)/(1 + 2*r);
	prob[1] = prob[2] = (1 - 4 * prob[0]) / 12;
}
template<> void genotypeProbabilitiesWithIntercross<8>(double (&prob)[3], int nAIGenerations, double r)
{
	double tmp = pow(1-r, nAIGenerations-1);
	prob[0] = (tmp *(1-r)*(1-r)*(1-r)/8 + (2*r + 1 - tmp)/64)/(1 + 2*r);
	prob[1] = prob[2] = (1 - 8 * prob[0]) / 56;
}
template<int nFounders> void genotypeProbabilitiesNoIntercross(double (&expandedProbabilities)[nFounders][nFounders], double r)
{
	double probabilities[3];
	genotypeProbabilitiesNoIntercross<nFounders>(probabilities, r);
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = 0; j < nFounders; j++)
		{
			expandedProbabilities[i][j] = probabilities[mask[i][j]];
		}
	}
}
template<int nFounders> void genotypeProbabilitiesWithIntercross(double (&expandedProbabilities)[nFounders][nFounders], int nAIGenerations, double r)
{
	double probabilities[3];
	genotypeProbabilitiesWithIntercross<nFounders>(probabilities, nAIGenerations, r);
	for(int i = 0; i < nFounders; i++)
	{
		for(int j = 0; j < nFounders; j++)
		{
			expandedProbabilities[i][j] = probabilities[mask[i][j]];
		}
	}
}
template<int nFounders, int maxMarkerAlleles> bool rfhaps_cpu_internal(rfhaps_cpu_args& args)
{
	int nFinals = args.finalsMatrix.nrow(), nRecombLevels = args.recombinationVector.size();
	int nDifferentFunnels = args.funnelEncodings.size();
	int marker2RangeSize = args.marker2End - args.marker2Start;
	int maxStart = std::max(args.marker1Start, args.marker2Start), minEnd = std::min(args.marker1End, args.marker2End);
	bool* allowableMarkerPatternsPtr = args.allowableMarkerPatterns.get();
	std::vector<double>& lineWeights = args.lineWeights;
	Rcpp::List finalDimNames = args.finalsMatrix.attr("dimnames");
	Rcpp::CharacterVector finalNames = finalDimNames[0];
	Rcpp::IntegerVector pedigreeIDColumn = Rcpp::as<Rcpp::IntegerVector>(args.pedigreeMatrix("id"));
	std::vector<int> fid = Rcpp::as<std::vector<int> >(args.fidVector);

	int nMarkerPatternIDs = args.markerEncodings.size();
	int maxAIGenerations = *std::max_element(args.nIntercrossingGenerations.begin(), args.nIntercrossingGenerations.end());
	
	typedef std::pair<markerPatternID, markerPatternID> markerPair;
	
	typedef std::array<std::array<double, maxMarkerAlleles>, maxMarkerAlleles> arrayType;
	typedef sharedArray<arrayType> PerAIGenerationData;
	typedef sharedArray<PerAIGenerationData> PerRecombinationFractionData;
	typedef std::map<markerPair, sharedArray<PerRecombinationFractionData> > PerMarkerPairData;
	
	//data structure goes marker pair -> recombination level -> number of AI generations -> funnelID -> marker1 value -> marker2Value. Basically a huge array. 
	PerMarkerPairData computedContributions;
#ifdef USE_OPENMP
	#pragma omp parallel for schedule(static, 1)
#endif
	//This is a big chunk of code, but does NOT grow with problem size. 
	for(int firstPattern = 0; firstPattern < nMarkerPatternIDs; firstPattern++)
	{
		int firstMarkerEncoding = args.markerEncodings[firstPattern];
		int firstMarkerPattern[nFounders];
		for(int i = 0; i < nFounders; i++)
		{
			firstMarkerPattern[i] = ((firstMarkerEncoding & (7 << (3*i))) >> (3*i));
		}
		//marker alleles have been encoded so they're of the form [0, nAlleles), so can just look for max value
		int nFirstMarkerAlleles = *std::max_element(firstMarkerPattern, firstMarkerPattern+nFounders);
		for(int secondPattern = 0; secondPattern < nMarkerPatternIDs; secondPattern++)
		{
			int secondMarkerEncoding = args.markerEncodings[secondPattern];
			int secondMarkerPattern[nFounders];
			for(int i = 0; i < nFounders; i++)
			{
				secondMarkerPattern[i] = ((secondMarkerEncoding & (7 << (3*i))) >> (3*i));
			}
			int nSecondMarkerAlleles = *std::max_element(secondMarkerPattern, secondMarkerPattern+nFounders);			
			
			markerPair currentPair(firstPattern, secondPattern);
			sharedArray<PerRecombinationFractionData> perRecombData(new PerRecombinationFractionData[nRecombLevels]);
			for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
			{
				double recombFraction = args.recombinationVector[recombCounter];
				PerRecombinationFractionData& currentRecomb = perRecombData[recombCounter];
				currentRecomb = PerRecombinationFractionData(new PerAIGenerationData[maxAIGenerations+1]);
				for(int intercrossingGeneration = 0; intercrossingGeneration <= maxAIGenerations; intercrossingGeneration++)
				{
					PerAIGenerationData& markerProbabilitiesAllFunnels = currentRecomb[intercrossingGeneration];
					markerProbabilitiesAllFunnels = PerAIGenerationData(new arrayType[nDifferentFunnels]);
					double haplotypeProbabilities[nFounders][nFounders];
					if(intercrossingGeneration == 0)
					{
						genotypeProbabilitiesNoIntercross<nFounders>(haplotypeProbabilities, recombFraction);
					}
					else
					{
						genotypeProbabilitiesWithIntercross<nFounders>(haplotypeProbabilities, intercrossingGeneration, recombFraction);
					}
					for(int funnelCounter = 0; funnelCounter < nDifferentFunnels; funnelCounter++)
					{
						arrayType& markerProbabilitiesThisFunnel = markerProbabilitiesAllFunnels[funnelCounter];
						memset(&(markerProbabilitiesThisFunnel[0][0]), 0, sizeof(double)*maxMarkerAlleles*maxMarkerAlleles);
						funnelEncoding enc = args.funnelEncodings[funnelCounter];
						int funnel[8];
						for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
						{
							funnel[founderCounter] = ((enc & (7 << (3*founderCounter))) >> (3*founderCounter));
						}
						
						for(int firstMarkerValue = 0; firstMarkerValue <= nFirstMarkerAlleles; firstMarkerValue++)
						{
							//firstFounder is the index of the founder within the current funnel
							for(int firstFounder = 0; firstFounder < nFounders; firstFounder++)
							{
								if(firstMarkerPattern[funnel[firstFounder]] == firstMarkerValue)
								{
									for(int secondMarkerValue = 0; secondMarkerValue <= nSecondMarkerAlleles; secondMarkerValue++)
									{
										for(int secondFounder = 0; secondFounder < nFounders; secondFounder++)
										{
											 if(secondMarkerPattern[funnel[secondFounder]] == secondMarkerValue)
											 {
												markerProbabilitiesThisFunnel[firstMarkerValue][secondMarkerValue] += haplotypeProbabilities[firstFounder][secondFounder];
											 }
										 }
									}
								}
							}
						}
						//now take logs of every value in markerProbabilities
						for(int firstMarkerValue = 0; firstMarkerValue <= nFirstMarkerAlleles; firstMarkerValue++)
						{
							for(int secondMarkerValue = 0; secondMarkerValue <= nSecondMarkerAlleles; secondMarkerValue++)
							{
								markerProbabilitiesThisFunnel[firstMarkerValue][secondMarkerValue] = log10(markerProbabilitiesThisFunnel[firstMarkerValue][secondMarkerValue]);
							}
						}
					}
				}
			}
			#pragma omp critical
			{
				computedContributions.insert(make_pair(currentPair, perRecombData));
			}
		}
	}
#ifdef USE_OPENMP
	#pragma omp parallel for schedule(static, 1)
#endif
	//This set of loops DOES grow with problem size. 
	for(int markerCounter1 = args.marker1Start; markerCounter1 < args.marker1End; markerCounter1++)
	{
		int markerPatternID1 = args.markerPatternIDs[markerCounter1];
		for(int markerCounter2 = args.marker2Start; markerCounter2 < args.marker2End; markerCounter2++)
		{
			//For some bits in the lower triangle we have already calculated a corresponding bit in the upper triangle, so don't recalculate these. This means we have some values in args.result which are never set, but later code is aware of this
			if(markerCounter2 >= maxStart && markerCounter2 < minEnd && markerCounter1 >= maxStart && markerCounter1 < minEnd && markerCounter2 < markerCounter1) continue;
			int markerPatternID2 = args.markerPatternIDs[markerCounter2];
			if(allowableMarkerPatternsPtr[markerPatternID1*nMarkerPatternIDs + markerPatternID2])
			{
				typename PerMarkerPairData::iterator patternData = computedContributions.find(markerPair(markerPatternID1, markerPatternID2));
				if(patternData == computedContributions.end()) throw std::runtime_error("Internal error");
				sharedArray<PerRecombinationFractionData> markerPairData = patternData->second;
				for(int recombCounter = 0; recombCounter < nRecombLevels; recombCounter++)
				{
					PerRecombinationFractionData& perRecombLevelData = markerPairData[recombCounter];
					
					for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
					{
						int marker1Value = args.finalsMatrix(finalCounter, markerCounter1);
						int marker2Value = args.finalsMatrix(finalCounter, markerCounter2);
						if(marker1Value != NA_INTEGER && marker2Value != NA_INTEGER)
						{
							int intercrossingGenerations = args.nIntercrossingGenerations[finalCounter];
							funnelID currentLineFunnelID = args.funnelIDs[finalCounter];
							arrayType& perMarkerGenotypeValues = perRecombLevelData[intercrossingGenerations][currentLineFunnelID];
							double contribution = perMarkerGenotypeValues[marker1Value][marker2Value];
							args.result[(long)(markerCounter1 - args.marker1Start) *(long)nRecombLevels*(long)marker2RangeSize + (long)(markerCounter2-args.marker2Start) * (long)nRecombLevels + (long)recombCounter] += lineWeights[finalCounter] * contribution;
						}
					}
				}
			}
		}
	}
	return true;
}
//here we transfer maxMarkerAlleles over to the templated parameter section - This can make a BIG difference to memory usage if this is smaller, and it's going into a type so it has to be templated.
template<int nFounders> bool rfhaps_cpu_internal(rfhaps_cpu_args& args)
{
	switch(args.maxMarkerAlleles)
	{
		case 1:
			return rfhaps_cpu_internal<nFounders, 1>(args);
		case 2:
			return rfhaps_cpu_internal<nFounders, 2>(args);
		case 3:
			return rfhaps_cpu_internal<nFounders, 3>(args);
		case 4:
			return rfhaps_cpu_internal<nFounders, 4>(args);
		case 5:
			return rfhaps_cpu_internal<nFounders, 5>(args);
		case 6:
			return rfhaps_cpu_internal<nFounders, 6>(args);
		case 7:
			return rfhaps_cpu_internal<nFounders, 7>(args);
		case 8:
			return rfhaps_cpu_internal<nFounders, 8>(args);
		default:
			throw std::runtime_error("Internal error");
	}
}
bool rfhaps_cpu(rfhaps_cpu_args& args)
{
	if(args.nFounders == 2)
	{
		return rfhaps_cpu_internal<2>(args);
	}
	else if(args.nFounders == 4)
	{
		return rfhaps_cpu_internal<4>(args);
	}
	else if(args.nFounders == 8)
	{
		return rfhaps_cpu_internal<8>(args);
	}
	else
	{
		Rprintf("Number of founders must be 2, 4 or 8\n");
		return false;
	}
}
