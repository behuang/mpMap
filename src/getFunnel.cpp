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
	else if(nFounders >= 8)
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
		if(nFounders == 16)
		{
		  int mmmRow = findIDInPedigree(mmmID, pedigree);
		  int mmfRow = findIDInPedigree(mmfID, pedigree);
		  int mfmRow = findIDInPedigree(mfmID, pedigree);
		  int mffRow = findIDInPedigree(mffID, pedigree);
		  int fmmRow = findIDInPedigree(fmmID, pedigree);
		  int fmfRow = findIDInPedigree(fmfID, pedigree);
		  int ffmRow = findIDInPedigree(ffmID, pedigree);
		  int fffRow = findIDInPedigree(fffID, pedigree);
		  int mmmmID = male(mmmRow), mmmfID=female(mmmRow), mmfmID=male(mmfRow), mmffID=female(mmfRow);
		  int mfmmID=male(mfmRow), mfmfID=female(mfmRow), mffmID=male(mffRow), mfffID=female(mffRow);
		  int fmmmID=male(fmmRow), fmmfID=female(fmmRow), fmfmID=male(fmfRow), fmffID=female(fmfRow);
		  int ffmmID=male(ffmRow), ffmfID=female(ffmRow), fffmID=male(fffRow), ffffID=female(fffRow);
		  funnel[0] = fidToFounderNumber(mmmmID);
		  funnel[1] = fidToFounderNumber(mmmfID);
		  funnel[2] = fidToFounderNumber(mmfmID);
		  funnel[3] = fidToFounderNumber(mmffID);
		  funnel[4] = fidToFounderNumber(mfmmID);
		  funnel[5] = fidToFounderNumber(mfmfID);
		  funnel[6] = fidToFounderNumber(mffmID);
		  funnel[7] = fidToFounderNumber(mfffID);
		  funnel[8] = fidToFounderNumber(fmmmID);
		  funnel[9] = fidToFounderNumber(fmmfID);
		  funnel[10] = fidToFounderNumber(fmfmID);
		  funnel[11] = fidToFounderNumber(fmffID);
		  funnel[12] = fidToFounderNumber(ffmmID);
		  funnel[13] = fidToFounderNumber(ffmfID);
		  funnel[14] = fidToFounderNumber(fffmID);
		  funnel[15] = fidToFounderNumber(ffffID);
		}
	}
	else
	{
		return false;
	}
	return true;
}
#undef fidToFounderNumber

