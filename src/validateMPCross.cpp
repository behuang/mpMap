#include "validateMPCross.h"
#include <Rinternals.h>
bool validateMPCross(Rcpp::RObject mpcross_, int& nFounders, std::string& error, bool checkPedigree, bool checkRF, bool checkLG, bool checkFID)
{
	Rcpp::Function asInteger("as.integer");
	Rcpp::IntegerVector dim;
	if(mpcross_.sexp_type() != VECSXP)
	{
		error = "Input mpcross object must have type list";
		return false;
	}
	Rcpp::List mpcross = Rcpp::as<Rcpp::List>(mpcross_);
	Rcpp::RObject classAttr_ = mpcross.attr("class");
	if(classAttr_.sexp_type() != STRSXP)
	{
		error = "Object class must be a character string";
		return false;
	}
	Rcpp::CharacterVector classAttr = mpcross.attr("class");
	if(classAttr.length() == 0)
	{
		error = "Internal error, class attribute had zero length";
		return false;
	}
	if(std::find(classAttr.begin(), classAttr.end(), "mpcross") == classAttr.end())
	{
		error = "Input objects must have class mpcross";
		return false;
	}
	Rcpp::CharacterVector mpcrossNames = mpcross.names();
	//check for existence of founders, finals, id
	if(std::find(mpcrossNames.begin(), mpcrossNames.end(), "founders") == mpcrossNames.end())
	{
		error = "Could not find entry mpcross$founders";
		return false;
	}
	if(std::find(mpcrossNames.begin(), mpcrossNames.end(), "finals") == mpcrossNames.end())
	{
		error = "Could not find entry mpcross$finals";
		return false;
	}
	if(std::find(mpcrossNames.begin(), mpcrossNames.end(), "id") == mpcrossNames.end())
	{
		error = "Could not find entry mpcross$id";
		return false;
	}
	if(checkPedigree)
	{
		if(std::find(mpcrossNames.begin(), mpcrossNames.end(), "pedigree") == mpcrossNames.end())
		{
			error = "Could not find entry mpcross$pedigree";
			return false;
		}
		Rcpp::RObject pedigree_ = mpcross["pedigree"];
		Rcpp::DataFrame pedigree;
		//Here we change the object in-place, if it's numeric instead of integer
		if(pedigree_.sexp_type() == REALSXP)
		{
			Rcpp::IntegerMatrix converted(Rf_coerceVector(pedigree_.get__(), INTSXP));
			DUPLICATE_ATTRIB(converted.get__(), pedigree_.get__());
			pedigree_ = mpcross["pedigree"] = converted;
		}
		//If it's a matrix, convert to data.frame
		if(pedigree_.sexp_type() == INTSXP)
		{
			//check it's 2-d
			if(!pedigree_.hasAttribute("dim"))
			{
				error = "mpcross$pedigree did not have a dimension attribute";
				return false;
			}
			dim = pedigree_.attr("dim");
			if(dim.length() != 2) 
			{
				error = "Internal error, pedigree had more than two dimensions";
				return false;
			}
			mpcross["pedigree"] = pedigree = Rcpp::Language("as.data.frame", pedigree_).eval();
		}
		else if(pedigree_.sexp_type() == VECSXP && Rcpp::as<std::string>(pedigree_.attr("class")) == "data.frame")
		{
			pedigree = Rcpp::as<Rcpp::DataFrame>(pedigree_);
		}
		else
		{
			error = "Input mpcross$pedigree must be an integer matrix or data.frame";
			return false;
		}
		//Colnames should be id, Male, Female, Observed, and then optionally Design
		Rcpp::List dimnames = Rcpp::Language("dimnames", pedigree).eval();
		Rcpp::IntegerVector dim = Rcpp::Language("dim", pedigree).eval();
		std::vector<std::string> colnames = Rcpp::as<std::vector<std::string> >(dimnames[1]);
		if(dim[1] != 4 && dim[1] != 5)
		{
			error = "Wrong number of columns in pedigree";
			return false;
		}
		if(colnames[0] != "id" || colnames[1] != "Male" || colnames[2] != "Female" || colnames[3] != "Observed")
		{
			error = "Column names of pedigree must be 'id', 'Male', 'Female' and 'Observed'";
			return false;
		}
		if(dim[1] == 5)
		{
			if(colnames[4] != "Design")
			{
                error = "If a pedigree has five columns the fifth must be named 'Design'";
				return false;
			}
		}
#define check_ped_column(INDEX, NAME)	\
		if(Rcpp::as<Rcpp::RObject>(pedigree(INDEX)).sexp_type() == REALSXP)	\
		{	\
			pedigree(INDEX) = asInteger(pedigree(INDEX));	\
		}	\
		if(Rcpp::as<Rcpp::RObject>(pedigree(INDEX)).sexp_type() != INTSXP)	\
		{	\
			error = "Pedigree column " NAME " must be numeric";	\
			return false;	\
		}	
		//convert first column, if it's numeric
		check_ped_column(0, "id");
		check_ped_column(1, "Male");
		check_ped_column(2, "Female");
		check_ped_column(3, "Observed");
#undef check_ped_column
		if(pedigree.length() == 5 && Rcpp::as<Rcpp::RObject>(pedigree("Design")).sexp_type() != STRSXP)
		{
			error = "Pedigree column Design must be a string";
			return false;
		}
	}
	//check type of founders
	Rcpp::RObject founders_ = mpcross["founders"];
	if(founders_.sexp_type() == REALSXP)
	{
		Rcpp::IntegerMatrix converted(Rf_coerceVector(founders_.get__(), INTSXP));
		DUPLICATE_ATTRIB(converted.get__(), founders_.get__());
		founders_ = mpcross["founders"] = converted;
	}
	if(founders_.sexp_type() != INTSXP)
	{
		error = "Input mpcross$founders must be an integer matrix";
		return false;
	}
	//check it's 2-d
	dim = founders_.attr("dim");
	if(dim.length() != 2) 
	{
		error = "Internal error, mpcross$founders had more than two dimensions";
		return false;
	}

	//check type of finals
	Rcpp::RObject finals_ = mpcross["finals"];
	if(finals_.sexp_type() == REALSXP)
	{
		Rcpp::IntegerMatrix converted(Rf_coerceVector(finals_.get__(), INTSXP));
		DUPLICATE_ATTRIB(converted.get__(), finals_.get__());
		finals_ = mpcross["finals"] = converted;
	}
	if(finals_.sexp_type() != INTSXP)
	{
		error = "Input mpcross$finals must be an integer matrix";
		return false;
	}
	//check it's 2-d
	dim = finals_.attr("dim");
	if(dim.length() != 2) 
	{
		error = "Internal error, mpcross$finals had more than two dimensions";
		return false;
	}
	Rcpp::IntegerMatrix finals = mpcross["finals"];
	Rcpp::IntegerMatrix founders = mpcross["founders"];
	
	//get final dim names
	Rcpp::List finalDimNames = finals_.attr("dimnames");
	//get founder dim names
	Rcpp::List founderDimNames = founders.attr("dimnames");
	//get number of rows
	nFounders = founders.nrow();
	int nMarkers = founders.ncol();
	//It's also acceptable to have NO rows in the founders / finals matrices, which happens when combining multiple densigns. 
	bool noGeneticData = Rcpp::as<Rcpp::RObject>(founderDimNames[0]).sexp_type() == NILSXP && Rcpp::as<Rcpp::RObject>(finalDimNames[0]).sexp_type() == NILSXP;
	//check that the column names are strings, and the row names are strings or there are no rows
	if(Rcpp::as<Rcpp::RObject>(founderDimNames[0]).sexp_type() != STRSXP && !noGeneticData)
	{
		error = "Row names of mpcross$founders must be character strings";
		return false;
	}
	if(Rcpp::as<Rcpp::RObject>(founderDimNames[1]).sexp_type() != STRSXP)
	{
		error = "Column names of mpcross$founders must be character strings";
		return false;
	}
	//check that dimension names are strings
	if(Rcpp::as<Rcpp::RObject>(finalDimNames[1]).sexp_type() != STRSXP)
	{
		error = "Column names of mpcross$finals must be character strings";
		return false;
	}
	if(Rcpp::as<Rcpp::RObject>(finalDimNames[0]).sexp_type() != STRSXP && !noGeneticData)
	{
		error = "Row names of mpcross$finals must be character strings";
		return false;
	}
	//check that founders and finals have same number of columns
	if(finals.ncol() != nMarkers)
	{
		error = "mpcross$finals and mpcross$founders had a different number of columns";
		return false;
	}
	std::vector<std::string> founderColNames = Rcpp::as<std::vector<std::string> >(founderDimNames[1]);
	std::vector<std::string> finalColNames = Rcpp::as<std::vector<std::string> >(finalDimNames[1]);
	if(!std::equal(founderColNames.begin(), founderColNames.end(), finalColNames.begin()))
	{
		error = "Column names of mpcross$founders and mpcross$finals must agree";
		return false;
	}
	if(checkFID)
	{
		if(std::find(mpcrossNames.begin(), mpcrossNames.end(), "fid") == mpcrossNames.end())
		{
			error = "Could not find entry mpcross$fid";
			return false;
		}
		//check type of fid
		Rcpp::RObject fid_ = mpcross["fid"];
		if(fid_.sexp_type() != INTSXP)
		{
			error = "Input mpcross$fid must be an integer vector";
			return false;
		}
		Rcpp::IntegerVector fid = mpcross["fid"];
		//confirm that this is same as number of entries as rows in founders
		if(nFounders != fid.length())
		{
			error = "Wrong number of entries in mpcross$fid";
			return false;
		}
	}
	if(checkRF)
	{
			if(std::find(mpcrossNames.begin(), mpcrossNames.end(), "rf") == mpcrossNames.end())
			{
				error = "Could not find entry mpcross$rf";
				return false;
			}
			Rcpp::RObject rf_ = mpcross["rf"];
			if(rf_.sexp_type() != VECSXP)
			{
				error = "Input mpcross$rf must be a list";
				return false;
			}
			Rcpp::List rf(static_cast<SEXP>(mpcross["rf"]));
			Rcpp::CharacterVector rfNames = rf.names();
			
			if(std::find(rfNames.begin(), rfNames.end(), "lod") == rfNames.end())
			{
				error = "Input mpcross must contain an entry named 'rf$lod'";
				return false;
			}
			
			//test lkhd, lod and theta
			Rcpp::RObject lod_ = rf["lod"];
			if(lod_.sexp_type() != REALSXP)
			{
				error = "Input mpcross$rf$lod must be numeric";
				return false;
			}
			if(lod_.hasAttribute("dim"))
			{
				Rcpp::IntegerVector lodDim = lod_.attr("dim");
				if(lodDim.length() != 2) 
				{
					error = "Input mpcross$rf$lod must be a matrix";
					return false;
				}
			}
			else
			{
				error = "Input mpcross$rf$lod must have a dimensions attribute";
				return false;
			}
			Rcpp::NumericMatrix lod = rf["lod"];
			if(lod.nrow() != lod.ncol() || lod.nrow() != nMarkers) 
			{
				error = "Input mpcross$rf$lod must be a square matrix";
				return false;
			}
			
			if(std::find(rfNames.begin(), rfNames.end(), "lkhd") == rfNames.end())
			{
				error = "Input mpcross must contain an entry named 'rf$lkhd'\n";
				return false;
			}
			
			Rcpp::RObject lkhd_ = rf["lkhd"];
			if(lkhd_.sexp_type() != REALSXP)
			{
				error = "Input mpcross$rf$lkhd must be numeric";
				return false;
			}
			if(lkhd_.hasAttribute("dim"))
			{
				Rcpp::IntegerVector lkhdDim = lkhd_.attr("dim");
				if(lkhdDim.length() != 2) 
				{
					error = "Input mpcross$rf$lkhd must be a matrix";
					return false;
				}
			}
			else
			{
				error = "Input mpcross$rf$lkhd must have a dimensions attribute";
				return false;
			}
			Rcpp::NumericMatrix lkhd = rf["lkhd"];
			if(lkhd.nrow() != lkhd.ncol() || lkhd.nrow() != nMarkers) 
			{
				error = "Input mpcross$rf$lkhd must be a square matrix";
				return false;
			}
			
			if(std::find(rfNames.begin(), rfNames.end(), "theta") == rfNames.end())
			{
				error = "Input mpcross must contain an entry named 'rf$theta'";
				return false;
			}
			
			Rcpp::RObject theta_ = rf["theta"];
			if(theta_.sexp_type() != REALSXP)
			{
				error = "Input mpcross$rf$theta must be numeric";
				return false;
			}
			if(theta_.hasAttribute("dim"))
			{
				Rcpp::IntegerVector thetaDim = theta_.attr("dim");
				if(thetaDim.length() != 2) 
				{
					error = "Input mpcross$rf$theta must be a matrix";
					return false;
				}
			}
			Rcpp::NumericMatrix theta = rf["theta"];
			if(theta.nrow() != theta.ncol() || theta.nrow() != nMarkers) 
			{
				error = "Input mpcross$rf$theta must be a square matrix";
				return false;
			}
			//check dimnames of all three matrices
			Rcpp::List thetaDimNames = theta.attr("dimnames");
			Rcpp::List lkhdDimNames = lkhd.attr("dimnames");
			Rcpp::List lodDimNames = lod.attr("dimnames");
			if(thetaDimNames.length() != 2 || lkhdDimNames.length() != 2 || lodDimNames.length() != 2)
			{
				error = "Internal error, objects had wrong number of dimensions";
				return false;
			}
			if(Rcpp::as<Rcpp::RObject>(thetaDimNames[1]).sexp_type() != STRSXP)
			{
				error = "Column names of mpcross$rf$theta must be a character vector";
				return false;
			}
			if(Rcpp::as<Rcpp::RObject>(thetaDimNames[0]).sexp_type() != STRSXP)
			{
				error = "Row names of mpcross$rf$theta must be a character vector";
				return false;
			}
			if(Rcpp::as<Rcpp::RObject>(lodDimNames[1]).sexp_type() != STRSXP)
			{
				error = "Column names of mpcross$rf$lod must be a character vector";
				return false;
			}
			if(Rcpp::as<Rcpp::RObject>(lodDimNames[0]).sexp_type() != STRSXP)
			{
				error = "Row names of mpcross$rf$lod must be a character vector";
				return false;
			}
			if(Rcpp::as<Rcpp::RObject>(lkhdDimNames[1]).sexp_type() != STRSXP)
			{
				error = "Column names of mpcross$rf$lkhd must be a character vector";
				return false;
			}
			if(Rcpp::as<Rcpp::RObject>(lkhdDimNames[0]).sexp_type() != STRSXP)
			{
				error = "Row names of mpcross$rf$lkhd must be a character vector";
				return false;
			}
			std::vector<std::string> thetaColumnNames = Rcpp::as<std::vector<std::string> >(thetaDimNames[1]);
			std::vector<std::string> lodColumnNames = Rcpp::as<std::vector<std::string> >(lodDimNames[1]);
			std::vector<std::string> lkhdColumnNames = Rcpp::as<std::vector<std::string> >(lkhdDimNames[1]);
			
			std::vector<std::string> thetaRowNames = Rcpp::as<std::vector<std::string> >(thetaDimNames[0]);
			std::vector<std::string> lodRowNames = Rcpp::as<std::vector<std::string> >(lodDimNames[0]);
			std::vector<std::string> lkhdRowNames = Rcpp::as<std::vector<std::string> >(lkhdDimNames[0]);
			if(!std::equal(founderColNames.begin(), founderColNames.end(), thetaRowNames.begin()))
			{
				error = "Row names of mpcross$rf$theta must match those of mpcross$founders";
				return false;
			}
			if(!std::equal(founderColNames.begin(), founderColNames.end(), lodRowNames.begin()))
			{
				error = "Row names of mpcross$rf$lod must match those of mpcross$founders";
				return false;
			}
			if(!std::equal(founderColNames.begin(), founderColNames.end(), lkhdColumnNames.begin()))
			{
				error = "Row names of mpcross$rf$lkhd must match those of mpcross$founders";
				return false;
			}
			
			if(!std::equal(founderColNames.begin(), founderColNames.end(), thetaColumnNames.begin()))
			{
				error = "Column names of mpcross$rf$theta must match those of mpcross$founders";
				return false;
			}
			if(!std::equal(founderColNames.begin(), founderColNames.end(), lodColumnNames.begin()))
			{
				error = "Column names of mpcross$rf$lod must match those of mpcross$founders";
				return false;
			}
			if(!std::equal(founderColNames.begin(), founderColNames.end(), lkhdColumnNames.begin()))
			{
				error = "Column names of mpcross$rf$lkhd must match those of mpcross$founders";
				return false;
			}
	}
	if(checkLG)
	{
		if(std::find(mpcrossNames.begin(), mpcrossNames.end(), "lg") == mpcrossNames.end())
		{
			error = "Could not find entry mpcross$lg";
			return false;
		}
		Rcpp::RObject lg_ = mpcross["lg"];
		if(lg_.sexp_type() != VECSXP)
		{
			error = "Input mpcross$lg must be a list";
			return false;
		}
		Rcpp::List lg(static_cast<SEXP>(mpcross["lg"]));
		Rcpp::CharacterVector lgNames = lg.names();
		if(std::find(lgNames.begin(), lgNames.end(), "all.groups") == lgNames.end())
		{
			error = "Input mpcross must contain an entry named 'lg$all.groups'";
			return false;
		}
		if(std::find(lgNames.begin(), lgNames.end(), "groups") == lgNames.end())
		{
			error = "Input mpcross must contain an entry named 'lg$groups'";
			return false;
		}
		//test all.groups and then groups
		Rcpp::RObject allGroups_ = lg["all.groups"];
		std::vector<int> allGroups;
		if(allGroups_.sexp_type() == REALSXP)
		{
			Rcpp::NumericVector tmp = allGroups_.get__();
			allGroups.resize(tmp.size());
			for(int i = 0; i < tmp.length(); i++)
			{
				allGroups[i] = tmp[i]+0.5f;
			}
		}
		else if (allGroups_.sexp_type() == INTSXP)
		{
			allGroups = Rcpp::as<std::vector<int> >(allGroups_);
		}
		else
		{
			error = "Input mpcross$lg$all.groups must be numeric";
			return false;
		}
		
		Rcpp::RObject groups_ = lg["groups"];
		std::vector<int> groups;
		if(groups_.sexp_type() == REALSXP)
		{
			Rcpp::NumericVector tmp = groups_.get__();
			groups.resize(tmp.size());
			for(int i = 0; i < tmp.length(); i++)
			{
				groups[i] = tmp[i]+0.5f;
			}
		}
		else if (groups_.sexp_type() == INTSXP)
		{
			groups = Rcpp::as<std::vector<int> >(groups_);
		}
		else
		{
			error = "Input mpcross$lg$groups must be numeric";
			return false;
		}
		if(groups.size() != (std::size_t)nMarkers)
		{
			error = "Input mpcross$lg$groups had the wrong length";
			return false;
		}
		
		Rcpp::RObject groupNames_ = Rcpp::as<Rcpp::RObject>(groups_.attr("names"));
		if(groupNames_.sexp_type() != STRSXP)
		{
			error = "Input mpcross$lg$groups must have a names attribute";
			return false;
		}
		std::vector<std::string> groupNames = Rcpp::as<std::vector<std::string> >(groupNames_);
		if(!std::equal(founderColNames.begin(), founderColNames.end(), groupNames.begin()))
		{
			error = "Names of mpcross$lg$groups must match column names of mpcross$founders";
			return false;
		}
		
		//groups must be contained in allGroups
		for(int i = 0; i < nMarkers; i++)
		{
			if(std::find(allGroups.begin(), allGroups.end(), groups[i]) == allGroups.end())
			{
				error = "All values of mpcross$lg$groups must be contained in mpcross$lg$all.groups";
				return false;
			}
		}
		//check that the groups vector contains contiguous chunks, having the same value. 
		for(std::size_t i = 0; i < allGroups.size(); i++)
		{
			int currentGroup = allGroups[i];
			
			std::vector<int>::iterator first = std::find(groups.begin(), groups.end(), currentGroup);
			if(first == groups.end()) continue;
			std::vector<int>::iterator last = std::find(groups.rbegin(), groups.rend(), currentGroup).base();
			while(first != last)
			{
				if(*first != currentGroup)
				{
					error = "All groups must appear as contiguous regions in mpcross$lg$groups. Please reorder mpcross using subset";
					return false;
				}
				first++;
			}
		}
	}
	return true;
}
RcppExport SEXP validate(SEXP object, SEXP checkPedigree_, SEXP checkRF_, SEXP checkLG_, SEXP checkFID_)
{
	Rcpp::RObject mpcross(object);
	bool checkPedigree = true, checkRF = true, checkLG = true, checkFID = true;
	try
	{
		Rcpp::IntegerVector checkPedigreeRcpp(checkPedigree_);
		checkPedigree = checkPedigreeRcpp[0];
	}
	catch(Rcpp::not_compatible&){}
	
	try
	{
		Rcpp::IntegerVector checkRFRcpp(checkRF_);
		checkRF = checkRFRcpp[0];
	}
	catch(Rcpp::not_compatible&){}
	
	try
	{
		Rcpp::IntegerVector checkLGRcpp(checkLG_);
		checkLG = checkLGRcpp[0];
	}
	catch(Rcpp::not_compatible&){}
	
	try
	{
		Rcpp::IntegerVector checkFIDRcpp(checkFID_);
		checkFID = checkFIDRcpp[0];
	}
	catch(Rcpp::not_compatible&){}

	int nFounders;
	std::string error;
	bool valid = validateMPCross(mpcross, nFounders, error, checkPedigree, checkRF, checkLG, checkFID);
	if(valid)
	{
		return R_NilValue;
	}
	return Rcpp::wrap(error);
}