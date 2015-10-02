#include "findIDInPedigree.h"
int findIDInPedigree(int id, Rcpp::DataFrame& pedigree)
{
	Rcpp::IntegerVector idVector = pedigree("id");
        if(id - 1 < pedigree.nrows() && idVector(id-1) == id)
        {
                return id-1;
        }
        else
        {
                int* i = std::find(&(idVector(0)), &(idVector(0)) + pedigree.nrows(), id);
                if(i == &(idVector(0)) + pedigree.nrows()) return -1;
                return i - &(idVector(0));
        }
}
