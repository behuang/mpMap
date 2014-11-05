#include "rfhaps.h"
#include <math.h>
#include "sharedArray.hpp"
#include <limits>
#include "validateMPCross.h"
#include <array>
#include "rfhapsSpecificDesign.h"
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
