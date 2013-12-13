#ifndef GET_FUNNEL_CPU_HEADER
#define GET_FUNNEL_CPU_HEADER
#include <Rcpp.h>
#include <vector>
bool getFunnel(int id, Rcpp::DataFrame pedigree, std::vector<int>& fid, int aiGenerations, int* funnel, int nPedigreeRows, int nFounders);
#endif
