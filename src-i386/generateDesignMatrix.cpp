#include "generateDesignMatrix.h"
SEXP generateDesignMatrix(SEXP n_, SEXP maxOffset_)
{
	char* stackmem;
	{
		std::string error;
		{
			Rcpp::RObject n_robject(n_);
			int n;
			if(n_robject.sexp_type() == REALSXP)
			{
				Rcpp::NumericVector tmp = n_;
				double tmp2 = Rcpp::as<double>(tmp);
				if(fabs(tmp2 - floor(tmp2)) > 1e-2) 
				{
					error = "Input n must be an integer";
					goto signal_error;
				}
				n = static_cast<int>(tmp2);
			}
			else if(n_robject.sexp_type() == INTSXP)
			{
				Rcpp::IntegerVector tmp = n_;
				n = Rcpp::as<int>(tmp);
			}
			else
			{
				error = "Input n must be an integer";
				goto signal_error;
			}

			Rcpp::RObject maxOffset_robject(maxOffset_);
			int maxOffset;
			if(maxOffset_robject.sexp_type() == REALSXP)
			{
				Rcpp::NumericVector tmp = maxOffset_;
				double tmp2 = Rcpp::as<double>(tmp);
				if(fabs(tmp2 - floor(tmp2)) > 1e-2) 
				{
					error = "Input maxOffset must be an integer";
					goto signal_error;
				}
				maxOffset = static_cast<int>(tmp2);
			}
			else if(maxOffset_robject.sexp_type() == INTSXP)
			{
				Rcpp::NumericVector tmp = maxOffset_;
				maxOffset = Rcpp::as<int>(tmp);
			}
			else
			{
				error = "Input n must be an integer";
				goto signal_error;
			}
			if(maxOffset > n)
			{
				error = "Input maxOffset cannot be larger than n";
				goto signal_error;
			}
			int resultRows = n*maxOffset - maxOffset*(maxOffset - 1)/2;
			Rcpp::IntegerMatrix result(resultRows, n);
			memset(&(result(0, 0)), 0, sizeof(int) * n * resultRows);
			
			for(int i = 0; i < n; i++)
			{
				int offset = 0;
				//j is the section going by rows
				for(int j = 0; j < maxOffset; j++)
				{
					//resultMat[offset + max(0, i-j):min(n-j, i-1) ,i] <- 1
					int end = std::min(n - j-1, i) + 1;
					for(int k = std::max(0, i-j); k < end; k++) result(offset + k, i) = 1;
					offset = offset + (n-j);
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