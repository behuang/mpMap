#include <stdio.h>
#include <stdlib.h>
#include <cuda.h>
#include <iostream>
#include <cassert>
#include "arrayDeleter.h"
#include "rfhaps_gpu.h"
#include <algorithm>
#include "R.h"
#include <exception>
#include "getFunnelGPU.h"
using namespace std;

#define SAFE_EXIT( m )\
  Rprintf("%s in file '%s' in line %i.\n", m,  __FILE__, __LINE__);\
  cudaDeviceReset();\
  exit(EXIT_FAILURE);

#define R_CUDA_SAFE_CALL( call )\
  {\
  cudaError_t cudaError = call ;\
  if( cudaError != cudaSuccess ) {\
    Rprintf("%s in file '%s' in line %i.\n", cudaGetErrorString(cudaError),  __FILE__, __LINE__);\
    cudaDeviceReset();\
    exit(EXIT_FAILURE);\
  }\
}

void selectGPU(int deviceNum) {
  int myDevice, numDevices;

  R_CUDA_SAFE_CALL( cudaGetDeviceCount( &numDevices ) );
  if (deviceNum > numDevices) {
    Rprintf("Unable to use device %i, only %i found.\n",deviceNum,numDevices);
  }
  if (deviceNum >= 0) {   
    // if caller specified a device then use it
    R_CUDA_SAFE_CALL( cudaSetDevice(deviceNum) );
  } else if (deviceNum == -1) {
    // take the first available (will share unless GPUs are in exclusive mode)
    Rprintf("Selecting first available GPU.\n");
    R_CUDA_SAFE_CALL( cudaSetDevice(0) );
  } else if (deviceNum == -2) {
    // try some smarts to round robin assign devices based on the MPI local rank
    char* cLocalRank;
    int localRank = 0;
    cLocalRank = getenv("OMPI_COMM_WORLD_LOCAL_RANK");
    if (cLocalRank!=NULL) {
      localRank = atoi(cLocalRank);
      Rprintf("Local rank is: %i.\n",localRank);
    } else {
      Rprintf("Unable to determine local rank.\n");
    }
    R_CUDA_SAFE_CALL( cudaSetDevice(localRank % numDevices) );
  } else {
    SAFE_EXIT("Unknown argument to selectGPU");
  }

  // which device did we end up with..
  R_CUDA_SAFE_CALL( cudaGetDevice( &myDevice ) );
  Rprintf("Using device %i.\n",myDevice);
}


template<int nFounders> __device__ void pr2pt(double r, double *prob);
template<> __device__ void pr2pt<4>(double r, double *prob)
{
	prob[0] = (1-r)/(4+8*r);
	prob[1] = r/(4+8*r);
	prob[2] = r/(4+8*r);
}
template<> __device__ void pr2pt<8>(double r, double *prob)
{
	prob[0] = (1-r)*(1-r)/(8+16*r);
	prob[1] = r*(1-r)/(8+16*r);
	prob[2] = r/(16+32*r);
}
template<int nFounders>__device__ void pr2ptirip(double r, int s, double *prob);
template<> __device__ void pr2ptirip<4>(double r, int s, double *prob)
{
  prob[0]=(pow(1-r, 2+s-1)/4+(2*r+1-pow(1-r, s-1))/16)/(1+2*r); 
  prob[1]=prob[2]=(1-4*prob[0])/12;
}
template<> __device__ void pr2ptirip<8>(double r, int s, double* prob)
{
        double tmp = pow(1-r, s-1);
        prob[0] = (tmp *(1-r)*(1-r)*(1-r)/8 + (2*r + 1 - tmp)/64)/(1 + 2*r);
        prob[1] = prob[2] = (1 - 8 * prob[0]) / 56;
}
extern __shared__ char dyn_shared_mem[]; /* dynamic allocated shared memory */

template<int nFounders>
__global__ void gpu_rfhaps(int nRecomb, int* ngen, 
			   int nPairs, int nFinals,
	   		   int *finalg, int* pair1, int* pair2, 
			   double *thvec, int* markerPatternIDs, bool* allowableMarkerPatterns, int nMarkerPatterns, double* lineWeights, double* output) {
        /*
	 *	Mask is a matrix that looks something like
	 *	   nfounders = 8      	nfounders = 4
	 *	   01222222		0122
 	 *	   10222222		1022
 	 *	   22012222		2201
 	 *	   22102222		2210
 	 *	   22220122
 	 *	   22221022
 	 *	   22222201
 	 *	   22222210
 	 */
	__shared__ int mask[8][8];
	int g1[8];
	int g2[8];
	/* I suppose this could be done by differently..
	 * by having each thread copy a portion of the mask from device memory
	 * seems pointless unless the mask gets very large
	 * I assume doing it this way just increases the code size
	 */
	mask[0][0]=0;	mask[0][1]=1;	mask[0][2]=2;	mask[0][3]=2;	mask[0][4]=2;	mask[0][5]=2;	mask[0][6]=2;	mask[0][7]=2;
	mask[1][0]=1;	mask[1][1]=0;	mask[1][2]=2;	mask[1][3]=2;	mask[1][4]=2;	mask[1][5]=2;	mask[1][6]=2;	mask[1][7]=2;
	mask[2][0]=2;	mask[2][1]=2;	mask[2][2]=0;	mask[2][3]=1;	mask[2][4]=2;	mask[2][5]=2;	mask[2][6]=2;	mask[2][7]=2;
	mask[3][0]=2;	mask[3][1]=2;	mask[3][2]=1;	mask[3][3]=0;	mask[3][4]=2;	mask[3][5]=2;	mask[3][6]=2;	mask[3][7]=2;
	mask[4][0]=2;	mask[4][1]=2;	mask[4][2]=2;	mask[4][3]=2;	mask[4][4]=0;	mask[4][5]=1;	mask[4][6]=2;	mask[4][7]=2;
	mask[5][0]=2;	mask[5][1]=2;	mask[5][2]=2;	mask[5][3]=2;	mask[5][4]=1;	mask[5][5]=0;	mask[5][6]=2;	mask[5][7]=2;
	mask[6][0]=2;	mask[6][1]=2;	mask[6][2]=2;	mask[6][3]=2;	mask[6][4]=2;	mask[6][5]=2;	mask[6][6]=0;	mask[6][7]=1;
	mask[7][0]=2;	mask[7][1]=2;	mask[7][2]=2;	mask[7][3]=2;	mask[7][4]=2;	mask[7][5]=2;	mask[7][6]=1;	mask[7][7]=0;
	double *shm_thvec = (double*)dyn_shared_mem; /* dynamically allocated shared memory */
	shm_thvec[threadIdx.x] = thvec[threadIdx.x];
	__syncthreads();

	/* work out which part of the r and k loops in the CPU implementation
	 * we are responsible for
	 */
	int r = blockIdx.x * blockDim.x + threadIdx.x;
	int k = blockIdx.y * blockDim.y + threadIdx.y;
	if (k>=nPairs) return;
	for (int i = 0; i < nFinals; i++) {
//	  assert(k >= 0 && i >= 0);
//	  assert(k < npairs && i < nfinals);
	  
	  /* TODO something seems to be corrupt in the finalg data
	   * need to find the cause of these rogue indiv values */

	  int p1 = pair1[k];
	  int p2 = pair2[k];
//	  assert(p1 >= 0 && p2 >= 0);
//	  assert(p1 < nmrk);
//	  assert(p2 < nmrk);

	  for (int j=0; j<nFounders; j++)
		g1[j] = g2[j] = 0;
	  
	  /* point to start of genotypes for the individual */
	  int h1 = finalg[p1*nFinals+i];
	  int h2 = finalg[p2*nFinals+i];
	  if ((h1>0)*(h2>0)) {  	/* check for missing values */
	    double theta = shm_thvec[r];
	    double probclass[3];
	    
	    if ((h1&1) == 1) {g1[0]=1; h1 -= 1; }
	    if ((h1&3) == 2) {g1[1]=1; h1 -= 2; }
	    if ((h1&7) == 4) {g1[2]=1; h1 -= 4; }
	    if ((h1&15) == 8){g1[3]=1; h1 -= 8; }
	    if ((h1&31) == 16){g1[4]=1; h1 -= 16;}
	    if ((h1&63) == 32){g1[5]=1; h1 -= 32;}
	    if ((h1&127) == 64){g1[6]=1; h1 -= 64;}
	    if ((h1&255) == 128){g1[7]=1;}

	    if ((h2&1) == 1) {g2[0]=1; h2 -= 1; }
	    if ((h2&3) == 2) {g2[1]=1; h2 -= 2; }
	    if ((h2&7) == 4) {g2[2]=1; h2 -= 4; }
	    if ((h2&15) == 8){g2[3]=1; h2 -= 8; }
	    if ((h2&31) == 16){g2[4]=1; h2 -= 16;}
	    if ((h2&63) == 32){g2[5]=1; h2 -= 32;}
	    if ((h2&127) == 64){g2[6]=1; h2 -= 64;}
	    if ((h2&255) == 128){g2[7]=1;}

	    /* Compute haplotype probabilities based on theta */
	    /* TODO template gpu_rfhaps based on ngen
	     * this will allow us to use 7 fewer registers
	     * for the case when ngen == 0
	     * and 2 fewer when ngen > 0
	     * it adds annoying complexity to the kernel invocation call portion though
	     */
		if(ngen[i] == 0)
		{
			pr2pt<nFounders>(theta, probclass);
		}
		else
		{
			pr2ptirip<nFounders>(theta, ngen[i], probclass);
		}
	    /* Check whether progeny genotypes are compatible with parent genotypes */
	    
	    /* For each combination of haplotypes which is feasible
	     * add the haplotype probabilities together */
	    double hp = 0;
	    for (int j=0; j<nFounders; j++){
	      for (int l=0; l<nFounders; l++){
		if (g1[j]*g2[l]) {
		  hp += probclass[mask[j][l]];
		}
	      }
	    }
	    /* log10(hp) is the individual contribution to the log-likelihood */
	    output[k*nRecomb+r] += (allowableMarkerPatterns[nMarkerPatterns * markerPatternIDs[p1] + markerPatternIDs[p2]] ? lineWeights[i]*log10(hp) : 0);
	    
	  } // end of check for missing values
  }
}
struct rfhaps_gpu_internal_args
{
	int* pair1, *pair2;
	long pairsOffset, nPairsToCalculate;
	int nMarkers, nFinals, recombOffset, nRecombToCalculate, nRecomb, nFounders;
	int* finalsD, *nIntercrossingD;
	double* recombinationFractionsD;
	bool hasAI;
	double* output;
	int marker2RangeSize, marker1Start, marker2Start;
	int* markerPatternIDs;
	bool* allowableMarkerPatterns;
	int nMarkerPatterns;
	double* lineWeightsD;
};
pedigreeColumns::pedigreeColumns(int* id, int* Male, int* Female, int* Observed, std::vector<std::string>& Design)
: id(id), Male(Male), Female(Female), Observed(Observed), Design(Design)
{}
__host__ bool rfhaps_gpu_internal(rfhaps_gpu_internal_args& args)
{
	assert(args.nFounders == 4 || args.nFounders == 8);
	const long resultSize = args.nPairsToCalculate * args.nRecombToCalculate;

	double* outputD;
	cudaMalloc((void**)&outputD, resultSize*sizeof(double));
	//cudaMemcpy(outputD, args.output, sizeof(double)*resultSize, cudaMemcpyHostToDevice);
	cudaMemset(outputD, 0, sizeof(double)*resultSize);


	//transfer pairs data
	int* pair1D, *pair2D;
	cudaMalloc((int**)&pair1D, args.nPairsToCalculate * sizeof(int));
	cudaMemcpy(pair1D, args.pair1 + args.pairsOffset, args.nPairsToCalculate*sizeof(int), cudaMemcpyHostToDevice);
	cudaMalloc((int**)&pair2D, args.nPairsToCalculate * sizeof(int));
	cudaMemcpy(pair2D, args.pair2 + args.pairsOffset, args.nPairsToCalculate*sizeof(int), cudaMemcpyHostToDevice);

	int threadsX = args.nRecombToCalculate;
	int threadsY = args.nPairsToCalculate;

	dim3 dimBlock(threadsX, floor(440/threadsX)); // logical max of 512 threads per block (only 440 per multiprocessor on Fermi anyway)
	dim3 dimGrid(1, ceil((double)threadsY / (double)dimBlock.y));
	size_t dynSharedSize = sizeof(double)*args.nRecombToCalculate;

	if (args.nFounders==4) 
	{
		gpu_rfhaps<4><<<dimGrid, dimBlock, dynSharedSize>>>(args.nRecombToCalculate, args.nIntercrossingD, args.nPairsToCalculate, args.nFinals, args.finalsD, pair1D, pair2D, args.recombinationFractionsD + args.recombOffset, args.markerPatternIDs, args.allowableMarkerPatterns, args.nMarkerPatterns, args.lineWeightsD, outputD);
	}
	else if (args.nFounders==8)
	{
		gpu_rfhaps<8><<<dimGrid, dimBlock, dynSharedSize>>>(args.nRecombToCalculate, args.nIntercrossingD, args.nPairsToCalculate, args.nFinals, args.finalsD, pair1D, pair2D, args.recombinationFractionsD + args.recombOffset, args.markerPatternIDs, args.allowableMarkerPatterns, args.nMarkerPatterns, args.lineWeightsD, outputD);
	}
	else
	{
		Rprintf("nFounders must have value 4 or 8\n");
		exit(-1);
	}
	bool result = true;
	cudaThreadSynchronize();
	cudaError_t lastError = cudaGetLastError();
	if(lastError != cudaSuccess)
	{
	  Rprintf("CUDA Last Error: %s\n",cudaGetErrorString(lastError));
		result = false;
	}

	double* copiedOutput = new double[resultSize];;
	cudaMemcpy(copiedOutput, outputD, resultSize*sizeof(double), cudaMemcpyDeviceToHost);

	for(int pairCounter = 0; pairCounter < args.nPairsToCalculate; pairCounter++)
	{
		int markerCounter2 = args.pair2[pairCounter + args.pairsOffset];
		int markerCounter1 = args.pair1[pairCounter + args.pairsOffset];
		for(int recombCounter = 0; recombCounter < args.nRecombToCalculate; recombCounter++)
		{
			//Turns out this overflows the range of a signed int
			long index = (long)(markerCounter1 - args.marker1Start)*(long)args.nRecomb*(long)args.marker2RangeSize + (long)(markerCounter2-args.marker2Start)*(long)args.nRecomb + (long)(recombCounter+args.recombOffset);
			args.output[index] += copiedOutput[(long)pairCounter*(long)args.nRecombToCalculate + (long)recombCounter];
		}
	}
	cudaFree(outputD);
	cudaFree(pair1D);
	cudaFree(pair2D);
	delete[] copiedOutput;
	return result;
}
extern "C" __host__ bool rfhaps_gpu(rfhaps_gpu_args& args)
{
        selectGPU(args.deviceNum);

	int marker2RangeSize = args.marker2End - args.marker2Start, marker1RangeSize = args.marker1End - args.marker1Start;
	int nMarkers = args.markerPatternIDs.size();
	int nMarkerPatterns = args.markerEncodings.size();
	const int finalsSize = nMarkers * args.nFinals;
	int* copiedFinals = new int[finalsSize];

	assert(args.nFounders == 4 || args.nFounders == 8);

	//working out the number of pairs is complicated because if we have a region on the diagonal we use the symmetry to avoid making double calculations. Whereas if we have a bit on the 
	//off-diagonal we need to calculate every value
	int maxStart = std::max(args.marker1Start, args.marker2Start);
	int minEnd = std::min(args.marker1End, args.marker2End);
	long square = std::max(minEnd - maxStart, 0);
	long squarePairs = square*(square + 1) /2;
	long nPairs = (marker1RangeSize * marker2RangeSize) - square * square + squarePairs;

	//re-encode finals genetic data, so that now the 1st bit says whether that individual is compatible with founder 1, 2nd bit compatible with founder 2, etc
	intArray8 funnel_;
	for(int individualCounter = 0; individualCounter < args.nFinals; individualCounter++)
	{
		funnel_ = args.funnels[individualCounter];
		for(int markerCounter = 0; markerCounter < nMarkers; markerCounter++)
		{
			int newValue = 0;
			int oldValue = args.finals[individualCounter+args.nFinals*markerCounter];
			for(int founderCounter = 0; founderCounter < args.nFounders; founderCounter++)
			{
				if(oldValue == args.founders[funnel_.val[founderCounter] - 1 + args.nFounders*markerCounter]) newValue += (1 << founderCounter);
			}
			copiedFinals[individualCounter+args.nFinals*markerCounter] = newValue;
		}
	}
	//transfer intercrossing data
	int* nIntercrossingD;
	cudaError_t cudaAllocResult = cudaMalloc((void**)&nIntercrossingD, args.nFinals * sizeof(int));
	if(cudaAllocResult != cudaSuccess)
	{
		Rprintf("Error calling cudaMalloc with %d bytes: %s\n", args.nFinals * sizeof(int), cudaGetErrorString(cudaAllocResult));
		return false;
	}
	else
	{
		Rprintf("Allocated %d bytes\n", args.nFinals * sizeof(int));
	}
	cudaMemcpy(nIntercrossingD, args.nIntercrossing, args.nFinals * sizeof(int), cudaMemcpyHostToDevice);

	int* pair1 = new int[nPairs], *pair2 = new int[nPairs];
	int* p1Ptr = pair1, *p2Ptr = pair2;
	//generate pairs
	for(int i = args.marker1Start; i < args.marker1End; i++)
	{
		for(int j = args.marker2Start; j < args.marker2End; j++)
		{
			if(i >= maxStart && i < minEnd && j >= maxStart && j < minEnd && j < i) continue;
			*p2Ptr = j;
			*p1Ptr = i;
			p1Ptr++; p2Ptr++;
		}
	}

	//copy across final genetic data
	int* finalsD;
	cudaAllocResult = cudaMalloc((void**)&finalsD, finalsSize * sizeof(int));
	if(cudaAllocResult != cudaSuccess)
	{
		Rprintf("Error calling cudaMalloc with %d bytes: %s\n", finalsSize * sizeof(int), cudaGetErrorString(cudaAllocResult));
		return false;
	}
	else
	{
		Rprintf("Allocated %d bytes\n", finalsSize * sizeof(int));
	}
	cudaMemcpy(finalsD, copiedFinals, finalsSize * sizeof(int), cudaMemcpyHostToDevice);

	delete[] copiedFinals;
	//copy across recombination fractions
	double* recombinationFractionsD;
	cudaAllocResult = cudaMalloc((void**)&recombinationFractionsD, args.nRecomb * sizeof(double));
	if(cudaAllocResult != cudaSuccess)
	{
		Rprintf("Error calling cudaMalloc with %d bytes: %s\n", args.nRecomb * sizeof(double), cudaGetErrorString(cudaAllocResult));
		return false;
	}
	else
	{
		Rprintf("Allocated %d bytes\n", args.nRecomb * sizeof(double));
	}

	cudaMemcpy(recombinationFractionsD, args.recombination, args.nRecomb * sizeof(double), cudaMemcpyHostToDevice);

	//copy across the allowable marker patterns data
	bool* allowableMarkerPatternsD;
	cudaAllocResult = cudaMalloc((void**)&allowableMarkerPatternsD, nMarkerPatterns * nMarkerPatterns * sizeof(bool));
	if(cudaAllocResult != cudaSuccess)
	{
		Rprintf("Error calling cudaMalloc with %d bytes: %s\n", nMarkerPatterns * nMarkerPatterns*sizeof(bool), cudaGetErrorString(cudaAllocResult));
		return false;
	}
	else
	{
		Rprintf("Allocated %d bytes\n", nMarkerPatterns * nMarkerPatterns * sizeof(bool));
	}
	cudaMemcpy(allowableMarkerPatternsD, args.allowableMarkerPatterns, nMarkerPatterns * nMarkerPatterns * sizeof(bool), cudaMemcpyHostToDevice);
	
	int* markerPatternIDsD;
	cudaMalloc((void**)&markerPatternIDsD, args.markerPatternIDs.size() * sizeof(int));
	cudaMemcpy(markerPatternIDsD, &(args.markerPatternIDs[0]), args.markerPatternIDs.size() * sizeof(int), cudaMemcpyHostToDevice);
	
	//copy across line weights data
	double* lineWeightsD;
	cudaMalloc((void**)&lineWeightsD, args.lineWeights.size() * sizeof(double));
	cudaMemcpy(lineWeightsD, &(args.lineWeights[0]), args.lineWeights.size() * sizeof(double), cudaMemcpyHostToDevice);
	
	int threadsX = args.nRecomb;
	int threadsY = nPairs;

	dim3 dimBlock(threadsX, floor(440/threadsX)); // logical max of 512 threads per block (only 440 per multiprocessor on Fermi anyway)
	dim3 dimGrid(1, ceil((double)threadsY / (double)dimBlock.y));

	int donePairs = 0;
	long pairsPerCall = nPairs;
	if(dimGrid.y > 65535)
	{
		pairsPerCall = 65535 * dimBlock.y;
		dimGrid.y = pairsPerCall/dimBlock.y;
		threadsY = pairsPerCall;
		Rprintf("Splitting into %ld cuda calls....\n", (long)((nPairs+pairsPerCall-1)/pairsPerCall));
	}
	int requiredThreads = threadsX*threadsY;
	int totalThreads = dimBlock.x*dimGrid.x*dimBlock.y*dimGrid.y;
	Rprintf("Total threads needed = %d\n",requiredThreads);
	Rprintf("Threads in grid = %d\n",totalThreads);
	Rprintf("Surplus threads = %d\n\n", totalThreads - requiredThreads); /* these will need to just sit idle */

	Rprintf("Threads per block  %d x %d = %d\n",dimBlock.x,dimBlock.y,dimBlock.x*dimBlock.y);
	Rprintf("Blocks in grid  %d x %d = %d\n",dimGrid.x,dimGrid.y,dimGrid.x*dimGrid.y);
	//END_DEBUG
	rfhaps_gpu_internal_args internal_args;
	internal_args.pair1 = pair1;
	internal_args.pair2 = pair2;
	internal_args.nMarkers = nMarkers;
	internal_args.nFinals = args.nFinals;
	internal_args.recombOffset = 0;
	internal_args.nRecombToCalculate = args.nRecomb;
	internal_args.nRecomb = args.nRecomb;
	internal_args.finalsD = finalsD;
	internal_args.nIntercrossingD = nIntercrossingD;
	internal_args.recombinationFractionsD = recombinationFractionsD;
	internal_args.hasAI = args.hasAI;
	internal_args.nFounders = args.nFounders;
	internal_args.output = args.output;
	internal_args.marker2RangeSize = marker2RangeSize;
	internal_args.marker1Start = args.marker1Start;
	internal_args.marker2Start = args.marker2Start;
	internal_args.markerPatternIDs = markerPatternIDsD;
	internal_args.allowableMarkerPatterns = allowableMarkerPatternsD;
	internal_args.nMarkerPatterns = nMarkerPatterns;
	internal_args.lineWeightsD = lineWeightsD;
	
	int counter = 0;
	while(donePairs < nPairs)
	{
		Rprintf("Making cuda call %d\n", counter+1);
		internal_args.pairsOffset = donePairs;
		if(donePairs + pairsPerCall >= nPairs)
		{
			internal_args.nPairsToCalculate = nPairs - donePairs;
		}
		else internal_args.nPairsToCalculate = pairsPerCall;
		bool result = rfhaps_gpu_internal(internal_args);
		if(!result)
		{
			Rprintf("A CUDA call failed, exiting...\n");
			return false;
		}
		donePairs += internal_args.nPairsToCalculate;
		counter++;
	}

	delete[] pair1;
	delete[] pair2;

	Rprintf("Finished all CUDA calls\n");
	cudaFree(finalsD);
	cudaFree(lineWeightsD);
	cudaFree(nIntercrossingD);
	cudaFree(recombinationFractionsD);
	return true;
}

