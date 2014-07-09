#include "getFunnelCPU.h"
#include "getAllFunnels.h"
#include "validateMPCross.h"
#include "getFunnelCPU.h"
#include "intercrossingGenerations.h"
#include <sstream>
SEXP getAllFunnels(SEXP Rmpcross)
{
	char* stackmem;
	{
		std::string error;
		{
			int nFounders;
			Rcpp::RObject mpcross_ = Rmpcross;
			bool valid = validateMPCross(mpcross_, nFounders, error, true, false, false);
			if(!valid)
			{
				goto signal_error;
			}
			Rcpp::List mpcross = Rmpcross;
			Rcpp::DataFrame pedigree(mpcross["pedigree"]);
			Rcpp::IntegerVector id = mpcross["id"];
			int nFinals = id.length();
			std::vector<int> fid = Rcpp::as<std::vector<int> >(mpcross["fid"]);
			Rcpp::IntegerMatrix output(id.length(), nFounders);
			std::vector<int> nIntercrossingGenerations;
			nIntercrossingGenerations.resize(nFinals, 0);
			//get number of intercrossing generations
			bool ok = intercrossingGenerations(pedigree, nFounders, id, nIntercrossingGenerations);
			if(!ok)
			{
				error = "Problem determining number of intercrossing generations";
				goto signal_error;
			}
			//now get the actual funnels from the pedigree
			int funnel[8];
			for(int i = 0; i < id.length(); i++)
			{
				ok = getFunnel(id[i], pedigree, fid, nIntercrossingGenerations[i], funnel, pedigree.nrows(), nFounders);
				if(!ok)
				{
					std::stringstream ss;
					ss << "Problem with pedigree, for individual number " << (i+1) << ", having id " << id[i];
					error = ss.str();
					goto signal_error;
				}
				for(int j = 0; j < nFounders; j++) output(i, j) = funnel[j];
			}
			return output;
		}
	signal_error:
		stackmem = (char*)alloca(error.size() + 4);
		memset(stackmem, 0, error.size() + 4);
		memcpy(stackmem, error.c_str(), error.size());
	}
	Rf_error(stackmem);
	return R_NilValue;
}
