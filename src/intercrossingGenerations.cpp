#include "intercrossingGenerations.h"
#include "findIDInPedigree.h"
bool intercrossingGenerations(Rcpp::DataFrame& pedigree, int nFounders, Rcpp::IntegerVector& mpcrossID, std::vector<int>& output)
{
	#define pedFind(id) findIDInPedigree(id, pedigree)
	int nFinals = mpcrossID.size();
	int nPedigreeRows = pedigree.nrows();
	Rcpp::IntegerVector male = Rcpp::as<Rcpp::IntegerVector>(pedigree("Male")), female = Rcpp::as<Rcpp::IntegerVector>(pedigree("Female"));
	for(int finalCounter = 0; finalCounter < nFinals; finalCounter++)
	{
		int currentPedRow = pedFind(mpcrossID[finalCounter]);
		if(currentPedRow < 0 || currentPedRow > nPedigreeRows) return false;

		//Counter to stop if the loop goes too long and might be infinite
		int loopCounter = 0;
		//Pick the last row and proceed backwards up the pedigree until we're through all the selfing generations
		while(male(currentPedRow) == female(currentPedRow))
		{
			int nextPedID = male(currentPedRow);
			if(nextPedID < 0 || nextPedID > nPedigreeRows) return false;
			currentPedRow = pedFind(nextPedID);

			if(currentPedRow < 0 || currentPedRow > nPedigreeRows) return false;

			loopCounter++;
			if(loopCounter > 2000) return false;
		}
		//When we reach an NA in the pedigree the while condition will terminate, which is an error
		if((male(currentPedRow) != male(currentPedRow)) || (female(currentPedRow) != female(currentPedRow)))
		{
			return false;
		}
		int ngen = 0;
		while(male(currentPedRow) > 0)
		{
			int nextPedID = male(currentPedRow);
			if(nextPedID < 0 || nextPedID > nPedigreeRows) return false;
			currentPedRow = pedFind(nextPedID);
			if(currentPedRow < 0 || currentPedRow > nPedigreeRows) return false;

			ngen++;
			if(ngen > 2000) return false;
		}
		//another check for NA values
		if(male(currentPedRow) != male(currentPedRow))
		{
			return false;
		}
		output[finalCounter] = ngen - (int)((log(nFounders) / log(2)) + 0.5);
		if(output[finalCounter] < 0) return false;
	}
	#undef pedFind
	return true;
}
