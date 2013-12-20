#include <Rcpp.h>
#include "findIDInPedigree.h"
#define fidToFounderNumber(id) std::distance(fid.begin(), std::find(fid.begin(), fid.end(), id)) + 1
bool getFunnel(int id, Rcpp::DataFrame pedigree, std::vector<int>& fid, int aiGenerations, int* funnel, int nPedigreeRows, int nFounders)
{
	Rcpp::IntegerVector idColumn = pedigree("id");
	int row = findIDInPedigree(id, pedigree);
	if(row == -1) return false;
	Rcpp::IntegerVector male = Rcpp::as<Rcpp::IntegerVector>(pedigree("Male")), female = Rcpp::as<Rcpp::IntegerVector>(pedigree("Female"));
	while(male(row) == female(row))

	{
		int newID = male(row);
		row = findIDInPedigree(newID, pedigree);
		if(row == -1) return false;
	}
	int counter = aiGenerations;
	while(counter)
	{
		int newID = female(row);
		row = findIDInPedigree(newID, pedigree);
		if(row == -1) return false; 
		counter--;
	}
	int motherID = female(row), fatherID = male(row);
	int motherRow = findIDInPedigree(motherID, pedigree), fatherRow = findIDInPedigree(fatherID, pedigree);
	int mmID = female(motherRow), mfID = male(motherRow), fmID = female(fatherRow), ffID = male(fatherRow);
	if(nFounders == 2)
	{
		funnel[0] = fidToFounderNumber(motherID);
		funnel[1] = fidToFounderNumber(fatherID);
	}
	else if(nFounders == 4)
	{
		funnel[0] = fidToFounderNumber(mmID);
		funnel[1] = fidToFounderNumber(mfID);
		funnel[2] = fidToFounderNumber(fmID);
		funnel[3] = fidToFounderNumber(ffID);
	}
	else if(nFounders == 8)
	{
		int mmRow = findIDInPedigree(mmID, pedigree), mfRow = findIDInPedigree(mfID, pedigree), fmRow = findIDInPedigree(fmID, pedigree), ffRow = findIDInPedigree(ffID, pedigree);
		int mmmID = female(mmRow), mmfID = male(mmRow), mfmID = female(mfRow), mffID = male(mfRow);
		int fmmID = female(fmRow), fmfID = male(fmRow), ffmID = female(ffRow), fffID = male(ffRow);
		funnel[0] = fidToFounderNumber(mmmID);
		funnel[1] = fidToFounderNumber(mmfID);
		funnel[2] = fidToFounderNumber(mfmID);
		funnel[3] = fidToFounderNumber(mffID);
		funnel[4] = fidToFounderNumber(fmmID);
		funnel[5] = fidToFounderNumber(fmfID);
		funnel[6] = fidToFounderNumber(ffmID);
		funnel[7] = fidToFounderNumber(fffID);
	}
	else
	{
		return false;
	}
	return true;
}
#undef fidToFounderNumber

