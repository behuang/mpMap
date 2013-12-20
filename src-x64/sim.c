#include "mpMap.h"

void creategam(double *genome, int nmrk, int *ibd, int *gamete)
{
  int i;
  int *chrm;
  double r;

  chrm = (int*) R_alloc(nmrk, sizeof(int));
  // do some setup
  r=Random(); 
  chrm[0] = 1;
  if (r<0.5) chrm[0] = 0;
  gamete[0] = ibd[chrm[0]*(nmrk)]; 
 
  for (i=1; i<nmrk; i++)
  {
    r=Random();
    if (r < genome[i-1]) chrm[i]=1-chrm[i-1];
    if (r > genome[i-1]) chrm[i]=chrm[i-1];
    gamete[i] = ibd[chrm[i]*(nmrk)+i];
  }

  return;
}

void gengeno(double *genome, double *genome2, int *id, int *mother, 
		int *father, int *nmrk, int *nfdr, int *nped, int *seed, 
		int *transpos, int *transval, int *out)
{
  int i, j;
  int *done, *ibdgen1, *ibdgen2, *gamete;
  int count=*nfdr;
  int prt1, prt2, prt1nm=0, prt2nm=0, sdi1, sdi2;

  done = (int*) R_alloc((*nped+*nfdr), sizeof(int));
  ibdgen1 = (int*) R_alloc(2*(*nmrk), sizeof(int));
  ibdgen2 = (int*) R_alloc(2*(*nmrk), sizeof(int));
  gamete = (int*) R_alloc(*nmrk, sizeof(int));

  // need to initialize the founders
  for (i=0; i<*nfdr; i++)
  for (j=0; j<2*(*nmrk); j++) {
 	out[i*2*(*nmrk)+j] = i+1;
	done[i] = i;
  }

  PlantSeeds(*seed); 

  for (i=*nfdr; i<*nped; i++)
  {
    // drop down through pedigree
    if (which(id[i]-1, done, count)<0)
    {
	prt1 = which(mother[i]-1, done, count);
	prt2 = which(father[i]-1, done, count);

	if (prt1>=0) prt1nm = done[prt1]+1;
	sdi1 = which(prt1nm, id, *nped);
	if (prt2>=0) prt2nm = done[prt2]+1;
	sdi2 = which(prt2nm, id, *nped);
	
	for (j=0; j<2*(*nmrk); j++)
	{
	  ibdgen1[j] = out[sdi1*2*(*nmrk)+j];
	  ibdgen2[j] = out[sdi2*2*(*nmrk)+j];
	}
	
	if (ibdgen1[*transpos]==*transval) 
		creategam(genome2, *nmrk, ibdgen1, gamete); 
	else 	creategam(genome, *nmrk, ibdgen1, gamete);
	for (j=0; j<(*nmrk); j++)
	  out[i*2*(*nmrk)+j] = gamete[j];

	if (ibdgen2[*transpos]==*transval)
		creategam(genome2, *nmrk, ibdgen2, gamete); 
	else	creategam(genome, *nmrk, ibdgen2, gamete);
	for (j=0; j<(*nmrk); j++)
	  out[i*2*(*nmrk)+j+(*nmrk)] = gamete[j];
	
	done[i]=i;
	count++;
     } // end if
  } // end for
  return; 
}
