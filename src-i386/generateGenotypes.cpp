#include "generateGenotypes.h"
#include <Rcpp.h>
void createGamete(Rcpp::NumericVector& recombinationFractions, Rcpp::IntegerVector& geneticData, int transFounder, Rcpp::IntegerVector& transMarkers, int* output)
{
	int nMarkers = recombinationFractions.length() + 1;
	if((transFounder == 0 && transMarkers.length() != 0) || (transFounder != 0 && transMarkers.length() == 0))
	{
		throw std::runtime_error("Internal error");
	}
	//if there's a translocation and it could be inherited through a gamete, slightly different code
	if(transFounder != 0) 
	{
		int minTransMarker = *std::min_element(transMarkers.begin(), transMarkers.end()), maxTransMarker = *std::max_element(transMarkers.begin(), transMarkers.end());
		if(maxTransMarker - minTransMarker != (transMarkers.length()-1) || maxTransMarker != transMarkers(transMarkers.length()-1) || minTransMarker != transMarkers(0)) throw std::runtime_error("Internal error");
		
		bool firstTrans = geneticData(transMarkers(0) - 1) == transFounder, secondTrans = geneticData(transMarkers(0) + nMarkers - 1) == transFounder;
		if(firstTrans || secondTrans)
		{
			//OK, we skew the fraction which inherit the translocation here. If both have the translocation then it makes no difference
			float cutOff;
			if(firstTrans && secondTrans) cutOff = 0.5;
			else if(firstTrans) cutOff = 7.0/12.0;
			else cutOff = 5.0/12.0;
			//float cutOff = 0.5;
			
			//if currentHaplotype == 0, we're using the first allele for the first translocation location
			//if currentHaplotype == 1 we're using the second allele
			int currentHaplotype = 0;
			if(Rcpp::as<float>(Rcpp::runif(1, 0, 1)) > cutOff) currentHaplotype = 1;
			//Put in the chunk of genetic data for the translocation
			for(Rcpp::IntegerVector::iterator i = transMarkers.begin(); i != transMarkers.end(); i++)
			{
				output[*i - 1] = geneticData[nMarkers * currentHaplotype + *i - 1];
			}
			//Put in the bits preceeding. We want to keep a copy of the haplotype at the translocation because we'll use it to go the other way too. 
			int runningHaplotype = currentHaplotype;
			for(int i = minTransMarker - 2; i != -1; i--)
			{
				//recombination if we fall below the relevant recombination fraction
				if(Rcpp::as<float>(Rcpp::runif(1, 0, 1)) < recombinationFractions(i)) runningHaplotype = abs(runningHaplotype-1);
				output[i] = geneticData[nMarkers * runningHaplotype + i];
			}
			//And then the bits following
			runningHaplotype = currentHaplotype;
			for(int i = maxTransMarker; i != nMarkers; i++)
			{
				if(Rcpp::as<float>(Rcpp::runif(1, 0, 1)) < recombinationFractions(i-1)) runningHaplotype = abs(runningHaplotype-1);
				output[i] = geneticData[nMarkers * runningHaplotype + i];
			}
			return;
		}
	}
	//Either there was no translocation originally specified, or it's not present in the individual we're producing gametes for
	//Set the first marker separately as it's set without reference to any other
	int runningHaplotype = 0;
	if(Rcpp::as<float>(Rcpp::runif(1, 0, 1)) < 0.5) runningHaplotype = 1;
	output[0] = geneticData[nMarkers * runningHaplotype];
	for(int i = 1; i < nMarkers; i++)
	{
		//recombination if we fall below the relevant recombination fraction
		if(Rcpp::as<float>(Rcpp::runif(1, 0, 1)) < recombinationFractions(i-1)) runningHaplotype = abs(runningHaplotype-1);
		output[i] = geneticData[nMarkers * runningHaplotype + i];
	}
}
SEXP generateGenotypes(SEXP RrecombinationFractions, SEXP Rpedigree, SEXP RtransMarkerIndices, SEXP RtransFounder)
{
	char* stackmem;
	{
		std::string error;
		{
			Rcpp::NumericVector recombinationFractions;
			Rcpp::DataFrame pedigree;
			Rcpp::IntegerVector transMarkerIndices;
			int transFounder;
			try
			{
				recombinationFractions = Rcpp::NumericVector(RrecombinationFractions);
			}
			catch(Rcpp::not_compatible&)
			{
				error = "Input recombinationDistances must be a numeric vector";
				goto signal_error;
			}
			try
			{
				pedigree = Rcpp::DataFrame(Rpedigree);
			}
			catch(Rcpp::not_compatible&)
			{
				error = "Input pedigree must be a data.frame";
				goto signal_error;
			}
			try
			{
				transMarkerIndices = Rcpp::IntegerVector(RtransMarkerIndices);
			}
			catch(Rcpp::not_compatible&)
			{
				error = "Input transMarkerIndices must be an integer vector";
				goto signal_error;
			}
			try
			{
				transFounder = Rcpp::as<int>(Rcpp::IntegerVector(RtransFounder));
			}
			catch(Rcpp::not_compatible&)
			{
				error = "Input transFounder must be a single integer value";
				goto signal_error;
			}
			//For no translocation we must have transFounder = 0 and no values in transMarkerIndices.
			//Otherwise we must have a non-zero value and have values in transMarkerIndices
			//anything else is an error
			if((transFounder == 0 && transMarkerIndices.length() != 0) || (transFounder != 0 && transMarkerIndices.length() == 0))
			{
				error = "Invalid combination of values for transFounder and transMarkerIndices. These must be either 0 and a vector of length 0, or a non-zero value and vector of non-zero length";
				goto signal_error;
			}
			//check transMarkerIndices
			if(transFounder != 0)
			{
				std::sort(transMarkerIndices.begin(), transMarkerIndices.end());
				for(int i = 0; i < transMarkerIndices.length()-1; i++)
				{
					if(transMarkerIndices(i) + 1 != transMarkerIndices(i+1))
					{
						error = "transMarkerIndices must be a set of contiguous markers";
						goto signal_error;
					}
				}
				//check valid indices
				for(int i = 0; i < transMarkerIndices.length(); i++)
				{
					if(transMarkerIndices(i) < 1 || transMarkerIndices(i) > recombinationFractions.size()+1) 
					{
						error = "translocation marker index was out of range";
						goto signal_error;
					}
				}
			}
			Rcpp::IntegerVector id = Rcpp::as<Rcpp::IntegerVector>(pedigree("id")), mother = Rcpp::as<Rcpp::IntegerVector>(pedigree("Female")), father = Rcpp::as<Rcpp::IntegerVector>(pedigree("Male"))	;
			int nMarkers = recombinationFractions.length() + 1, nPedRows = pedigree.nrows();
			//Columns 0 and nMarkers correspond to the the pair of alleles at the first marker
			Rcpp::IntegerMatrix result(pedigree.nrows(), 2*nMarkers);
			//once we encounter a line with parents which are not set to zero we set this to true. And subsequently don't allow any zero-values for mother or father
			bool finishedFounders = false;
			int founderCounter = 0;
			//count number of founders by looking at the number of rows with mother and founder both set to '0'. This DOESN'T have to be a power of 2. It just indicates the number of lines which are generated without reference to any parent lines.
			//So they're generated a little differently. 
			for(int lineCounter = 0; lineCounter < nPedRows; lineCounter++)
			{
				bool isFounder = (mother(lineCounter) == 0 || father(lineCounter) == 0);
				if(finishedFounders && (isFounder))
				{
					error = "Founder lines must be placed at the top of the pedigree";
					goto signal_error;
				}
				if(isFounder)
				{
					//If it's a founder, fill it with the index of the founder
					for(int counter = 0; counter < 2*nMarkers; counter++)
					{
						result(lineCounter, counter) = founderCounter+1;
					}
					founderCounter++;
				}
				else
				{
					//otherwise use the previously generated genetic data for the mother and father
					int motherPedigreeRowIndex = std::distance(id.begin(), std::find(id.begin(), id.end(), mother(lineCounter)));
					int fatherPedigreeRowIndex = std::distance(id.begin(), std::find(id.begin(), id.end(), father(lineCounter)));
					//copy out the genetic data for mother and father (includes BOTH alleles at every location)
					Rcpp::IntegerVector motherData = result.row(motherPedigreeRowIndex);
					Rcpp::IntegerVector fatherData = result.row(fatherPedigreeRowIndex);
					
					Rcpp::IntegerVector newGenotypes(2*nMarkers);
					createGamete(recombinationFractions, motherData, transFounder, transMarkerIndices, &(newGenotypes(0)));
					createGamete(recombinationFractions, fatherData, transFounder, transMarkerIndices, &(newGenotypes(0)) + nMarkers);
					result.row(lineCounter) = newGenotypes;
				}
			}
			return result;
		}
	signal_error:
		stackmem = (char*)alloca(error.size() + 4);
		memset(stackmem, 0, error.size() + 4);
		memcpy(stackmem, error.c_str(), error.size());
	}
	Rf_error(stackmem);
	return R_NilValue;
}