#ifndef _RFHAPS_GPU_H
#define _RFHAPS_GPU_H
#include <vector>
#include <string.h>
#include "Unique.hpp"
#include "rfhaps_common.h"
struct pedigreeColumns
{
	pedigreeColumns(int* id, int* Male, int* Female, int* Observed, std::vector<std::string>& Design);
	int* id;
	int* Male;
	int* Female;
	int* Observed;
	std::vector<std::string> Design;
};
struct rfhaps_gpu_args
{
	rfhaps_gpu_args(pedigreeColumns& pedigree, std::vector<markerPatternID>& markerPatternIDs, std::vector<int>& fidVector, std::vector<double>& lineWeights, std::vector<markerEncoding>& markerEncodings, std::vector<funnelID>& funnelIDs, std::vector<funnelEncoding>& funnelEncodings)
	: pedigree(pedigree), markerPatternIDs(markerPatternIDs), fid(fidVector), lineWeights(lineWeights), markerEncodings(markerEncodings), funnelIDs(funnelIDs), funnelEncodings(funnelEncodings)
	{}
	pedigreeColumns& pedigree;
	std::vector<markerPatternID>& markerPatternIDs;
	std::vector<int>& fid;
	int* founders;
	int* finals;
	int* fidVector;
	
	int nPedigreeRows;
	double* recombination;
	int* nIntercrossing;
	int nFounders;
	int nFinals;
	int* IDs;
	int nRecomb;
	bool hasAI;
	//this has already been zeroed by the caller and must be added too, NOT just overwritten. May be used again, in the case of multiple designs examined jointly. 
	double* output;
	int marker1Start, marker1End;
	int marker2Start, marker2End;
	bool* allowableMarkerPatterns;
	std::vector<double>& lineWeights;
	std::vector<markerEncoding>& markerEncodings;
	std::vector<funnelID>& funnelIDs;
	std::vector<funnelEncoding>& funnelEncodings;
        int deviceNum;
};
#ifdef HAS_CUDA
extern "C" bool rfhaps_gpu(rfhaps_gpu_args& args);
#endif
#endif
