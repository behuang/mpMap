#include "rfhaps.h"
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
#include <tr1/array>
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
bool rfhapsSpecificDesign(SEXP finals, SEXP founders, SEXP pedigree_, SEXP id, SEXP fid_, SEXP recombinationFractions, long marker1Start, long marker1End, long marker2Start, long marker2End, std::vector<double>& lineWeights, SEXP RuseGPU, SEXP RdeviceNum, std::tr1::shared_ptr<double> result, std::string* error)
{
	Rcpp::IntegerMatrix foundersMatrix(founders);
	Rcpp::IntegerMatrix finalsMatrix(finals);
	Rcpp::DataFrame pedigreeDataFrame(pedigree_);
	Rcpp::NumericVector recombinationVector(recombinationFractions);
	Rcpp::IntegerVector useGPU_(RuseGPU);
	Rcpp::IntegerVector deviceNum_(RdeviceNum);
	Rcpp::IntegerVector IDVector(id);
	Rcpp::IntegerVector fidVector(fid_);

	Rcpp::List finalsDimnames = finalsMatrix.attr("dimnames");
	SEXP markerNamesSEXP = finalsDimnames[1];
	if(TYPEOF(markerNamesSEXP) != STRSXP)
	{
		*error = "Input finals must have columns names";
		return false;
	}
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
	//Check the funnels here first, where it's easier to back out and report errors. We check it later too, but the cost is so minimal we can afford to do it twice 
	int messages = 0;
	std::vector<intArray8> allFunnels;
	for(long individualCounter = 0; individualCounter < nFinals; individualCounter++)
	{
		//if we only had 4 founders then there second half of this is garbage, we don't use that bit...
		intArray8 funnel, copiedFunnel;
		bool ok = getFunnel(IDVector(individualCounter), pedigreeDataFrame, fid, nIntercrossingGenerations[individualCounter], &(funnel.val[0]), pedigreeDataFrame.nrows(), nFounders);
		if(!ok)
		{
			Rprintf("Problem with pedigree for individual %d\n", individualCounter+1);
		}
		//check for individuals with funnels that don't contain all founders
		memcpy(&copiedFunnel, &funnel, sizeof(intArray8));
		std::sort(&(copiedFunnel.val[0]), &(copiedFunnel.val[0]) + nFounders);
		if(std::unique(&(copiedFunnel.val[0]), &(copiedFunnel.val[0]) + nFounders) != &(copiedFunnel.val[0]) + nFounders && messages < 6)
		{
			//This is only a warning by itself....
			Rprintf("Warning: Funnel for individual %d contained founders {%d, ", individualCounter+1, funnel.val[0]);
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
			//....but is more serious if it causes the observed marker data to be impossible. So check for this.
			for(int markerCounter = 0; markerCounter < nMarkers; markerCounter++)
			{
				bool okMarker = false;

				//If observed value is an NA then than's ok, continue
				int value = finalsMatrix(individualCounter, markerCounter);
				if(value == NA_INTEGER) continue;

				for(int founderCounter = 0; founderCounter < nFounders; founderCounter++)
				{
					if(finalsMatrix(individualCounter, markerCounter) == foundersMatrix(funnel.val[founderCounter]-1, markerCounter))
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
			if(messages == 6) Rprintf("Suppressing further output\n");
		}
		orderFunnel(&(funnel.val[0]), nFounders);
		allFunnels.push_back(funnel);
	}
	//re-code the founder and final marker genotypes so that they always start at 0 and go up to n-1 where n is the number of distinct marker alleles
	//We do this to make it easier to identify markers with identical segregation patterns. recodedFounders = column major matrix
	Rcpp::IntegerMatrix recodedFounders(nFounders, nMarkers), recodedFinals(nFinals, nMarkers);
	recodedFinals.attr("dimnames") = finalsMatrix.attr("dimnames");
	unsigned int maxAlleles = 0;
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
	if(maxAlleles > 8) 
	{
		*error = "Internal error - Cannot have more than eight alleles per marker";
		return false;
	}
	//map containing encodings of all the haplotypes (encoded into the key), and an associated unique index (we can't use the haplotype encoding as an index because they'll jump around a lot and could be quite big values sometimes I think). Unique indices are guaranteed to be contiguous numbers starting from 0 - So [0, markerPatterns.size()]. 
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
SEXP rfhaps(SEXP RmpcrossList, SEXP recombinationFractions, SEXP marker1Range, SEXP marker2Range, SEXP RlineWeights, SEXP RuseGPU, SEXP RdeviceNum)
{
	//Necessary because of the retarded longjmp use inside Rf_error
	char* stackmem;
	{
		std::string error;
		{
			Rcpp::RObject recombinationVector_(recombinationFractions);
			if(recombinationVector_.sexp_type() != REALSXP)
			{
				error = "Input recombinationFractions must be numeric";
				goto signal_error;
			}
			Rcpp::NumericVector recombinationVector(recombinationFractions);
			long nRecombLevels = recombinationVector.size();
			//must have 0.5 as one of the recombination fractions to consider
			int halfIndex = -1;
			for(int i = 0; i < nRecombLevels; i++)
			{
				if(recombinationVector(i) == 0.5) halfIndex = i;
			}
			if(halfIndex == -1) 
			{
				error = "Recombination vector must contain a value of 0.5";
				goto signal_error;
			}

			Rcpp::RObject mpcrossList_(RmpcrossList);
			if(mpcrossList_.sexp_type() != VECSXP)
			{
				error = "Internal Error: First input to .Call(\"rfhaps\", ...) must be a list";
				goto signal_error;
			}
			Rcpp::List mpcrossList(RmpcrossList);
			if(mpcrossList.length() == 0)
			{
				error = "Internal Error: First input to .Call(\"rfhaps\", ...) cannot have length 0";
				goto signal_error;
			}
			
			Rcpp::RObject lineWeights_(RlineWeights);
			if(lineWeights_.sexp_type() != VECSXP)
			{
				error = "Input lineWeights must be a list";
				goto signal_error;
			}
			Rcpp::List lineWeights(RlineWeights);
			if(lineWeights.length() != mpcrossList.length())
			{
				error = "Input lineWeights had the wrong length";
				goto signal_error;
			}
			std::vector<std::string> markerNames;
			for(int i = 0; i < mpcrossList.length(); i++)
			{
				int nFounders;
				bool valid = validateMPCross(mpcrossList[i], nFounders, error, true, false, false);
				if(!valid) goto signal_error;
				//Check line weights type
				Rcpp::RObject lineWeightsThisDesign_ = Rcpp::as<Rcpp::RObject>(lineWeights[i]);
				if(lineWeightsThisDesign_.sexp_type() != REALSXP)
				{
					error = "Input lineWeights contained object of the wrong type";
					goto signal_error;
				}
				std::vector<double> lineWeightsThisDesign = Rcpp::as<std::vector<double> >(lineWeightsThisDesign_);

				//Check marker names are the same for each input object
				Rcpp::List mpcross = mpcrossList[i];
				Rcpp::List dimnames = Rcpp::as<Rcpp::RObject>(mpcross["finals"]).attr("dimnames");
				if(dimnames.size() != 2)
				{
					error = "Finals matrix must have row and column names";
					goto signal_error;
				}
				Rcpp::CharacterVector finalColNames = dimnames[1];
				if(markerNames.size() == 0)
				{
					markerNames = Rcpp::as<std::vector<std::string> >(finalColNames);
				}
				else
				{
					//check for same number of markers and same names as in the first object
					std::vector<std::string> otherMarkerNames = Rcpp::as<std::vector<std::string> >(finalColNames);
					if(otherMarkerNames.size() != markerNames.size())
					{
						error = "Inconsistent marker names (colnames(mpcross$finals)) in the input mpcross objects";
						goto signal_error;
					}
					if(!std::equal(markerNames.begin(), markerNames.end(), otherMarkerNames.begin()))
					{
						error = "Inconsistent marker names (colnames(mpcross$finals)) in the input mpcross objects";
						goto signal_error;
					}
				}
				//check that the length of lineWeights matches up with the number of rows in the finals matrix
				unsigned int nFinals = Rcpp::as<Rcpp::IntegerMatrix>(mpcross["finals"]).nrow();
				if(lineWeightsThisDesign.size() != nFinals)
				{
					error = "An entry in lineWeights had the wrong length";
					goto signal_error;
				}
			}
			//check the marker range inputs
			long nMarkers = markerNames.size();
			int marker1Start, marker1End, marker2Start, marker2End;
			if(TYPEOF(marker1Range) == NILSXP)
			{
				marker1Start = 0;
				marker1End = nMarkers;
			}
			else
			{
				Rcpp::NumericVector range(marker1Range);
				if(range.size() != 2) 
				{
					error = "Input marker1Range must have length 2";
					goto signal_error;
				}
				marker1Start = range(0)-1;
				marker1End = range(1)-1;
			}

			if(TYPEOF(marker2Range) == NILSXP)
			{
				marker2Start = 0;
				marker2End = nMarkers;
			}
			else
			{
				Rcpp::NumericVector range(marker2Range);
				if(range.size() != 2) 
				{
					error = "Input marker2Range must have length 2";
					goto signal_error;
				}
				marker2Start = range(0)-1;
				marker2End = range(1)-1;
			}
			if(marker1Start >= marker1End || marker1Start < 0 || marker1End > nMarkers) 
			{
				error = "Invalid range passed to rfhaps for marker 1";
				goto signal_error;
			}
			if(marker2Start >= marker2End || marker2Start < 0 || marker2End > nMarkers) 
			{
				error = "Invalid range passed to rfhaps for marker 2";
				goto signal_error;
			}
			long marker1RangeSize = marker1End - marker1Start, marker2RangeSize = marker2End - marker2Start;
			//Indexing has form result[markerCounter1 *nRecombLevels*nMarkers + markerCounter2 * nRecombLevels + recombCounter]
			//NOT of the form result[recombCounter *nMarkers*nMarkers + markerCounter2 * nMarkers + markerCounter1], which is the form to be used if this was going to be processed in R
			//This is not an Rcpp::NumericVector because it can quite easily overflow the size of such a vector (signed int)
			sharedArray<double> result(new double[marker1RangeSize * marker2RangeSize * nRecombLevels]);
			memset((void*)result.get(), 0, marker1RangeSize * marker2RangeSize * nRecombLevels * sizeof(double));

			for(int i = 0; i < mpcrossList.length(); i++)
			{
				Rcpp::List mpcross = mpcrossList[i];
				std::vector<double> lineWeightsThisDesign = Rcpp::as<std::vector<double> >(lineWeights[i]);
				
				bool successful = rfhapsSpecificDesign(mpcross["finals"], mpcross["founders"], mpcross["pedigree"], mpcross["id"], mpcross["fid"], recombinationFractions, marker1Start, marker1End, marker2Start, marker2End, lineWeightsThisDesign, RuseGPU, RdeviceNum, result, &error);
				if(!successful) goto signal_error;
			}
			//now for some post-processing
			Rcpp::NumericMatrix theta(marker1RangeSize, marker2RangeSize), lod(marker1RangeSize, marker2RangeSize), lkhd(marker1RangeSize, marker2RangeSize);
			
			Rcpp::CharacterVector markerNames1(marker1RangeSize), markerNames2(marker2RangeSize);
			std::copy(markerNames.begin() + marker1Start, markerNames.begin() + marker1Start + marker1RangeSize, markerNames1.begin());
			std::copy(markerNames.begin() + marker2Start, markerNames.begin() + marker2Start + marker2RangeSize, markerNames2.begin());
			
			Rcpp::List outputDimNames(Rcpp::List::create(markerNames1, markerNames2));
			theta.attr("dimnames") = lod.attr("dimnames") = lkhd.attr("dimnames") = outputDimNames;

			double* resultPtr = result.get();
			int maxStart = std::max(marker1Start, marker2Start);
			int minEnd = std::min(marker1End, marker2End);
			//first get maximum and minimum likelihoods
#ifdef USE_OPENMP
			#pragma omp parallel for
#endif
			for(int markerCounter1 = marker1Start; markerCounter1 < marker1End; markerCounter1++)
			{
				for(int markerCounter2 = marker2Start; markerCounter2 < marker2End; markerCounter2++)
				{
					int destIndex1 = markerCounter1 - marker1Start, destIndex2 = markerCounter2 - marker2Start;
					int sourceIndex1 = destIndex1, sourceIndex2 = destIndex2;
					if(markerCounter2 >= maxStart && markerCounter2 < minEnd && markerCounter1 >= maxStart && markerCounter1 < minEnd && markerCounter2 < markerCounter1)
					{
						sourceIndex1 = markerCounter2 - maxStart + (maxStart - marker1Start);
						sourceIndex2 = markerCounter1 - maxStart + (maxStart - marker2Start);
					}

					double* start = resultPtr + (long)sourceIndex1 *(long)marker2RangeSize*(long)nRecombLevels + (long)sourceIndex2 * (long)nRecombLevels;
					double* end = resultPtr + (long)sourceIndex1 *(long)marker2RangeSize*(long)nRecombLevels + (long)sourceIndex2 * (long)nRecombLevels + (long)nRecombLevels;
					double* maxPtr = std::max_element(start, end), *minPtr = std::min_element(start, end);
					double max = *maxPtr, min = *minPtr;
					double currentTheta, currentLod;
					//This is the case where no data was available, across any of the experiments. This is precise, no numerical error involved
					if(max == 0 && min == 0)
					{
						currentLod = currentTheta = std::numeric_limits<double>::quiet_NaN();
					}
					else
					{
						currentTheta = recombinationVector(maxPtr - start);
						currentLod = max - resultPtr[(long)sourceIndex1*(long)nRecombLevels*(long)marker2RangeSize + (long)sourceIndex2 * (long)nRecombLevels + halfIndex];
					}
					theta(destIndex1, destIndex2) = currentTheta;
					lkhd(destIndex1, destIndex2) = max;
					lod(destIndex1, destIndex2) = currentLod;
				}
			}
			return Rcpp::List::create(Rcpp::Named("theta") = theta, Rcpp::Named("lod") = lod, Rcpp::Named("lkhd") = lkhd, Rcpp::Named("r") = recombinationVector);
		}
signal_error:
		stackmem = (char*)alloca(error.size() + 4);
		memset(stackmem, 0, error.size() + 4);
		memcpy(stackmem, error.c_str(), error.size());
	}
	Rf_error(stackmem);
	return R_NilValue;
}
