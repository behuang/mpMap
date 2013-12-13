#include "mpMap.h"

// set this up so that i can run a scan over a grid of positions, calculate
// likelihood at each point. 
void add3pt(int *finalg, int *founderg, int *id, int *mother, int *father, int *nfinals, int *nfounders, int *nmrk, int *nped, int* ngen, double *r12, double *r23, double *r13, int *left, int *mid, int *right, int *totmrk, int *npos, int *out, double *lod)
{
 int i, j, k, l, m, q, lk, rk, mk;
 int *funnel, *geni, *hapi;
 double *hprob, *p, *x, *hproblod, *plod, *xlod;
 int *cf1, *cf2, *cf3;
 int *fou, np=*npos/2;
 double hp;
 double *lkhd, *lkhd2;

 fou = (int*) R_alloc(*nfounders*3, sizeof(int));
 p = (double*) R_alloc(10, sizeof(double));
 plod = (double*) R_alloc(10, sizeof(double));
 x = (double*) R_alloc(10, sizeof(double));
 xlod = (double*) R_alloc(10, sizeof(double));
 hprob = (double*) R_alloc(20, sizeof(double));
 hproblod = (double*) R_alloc(20, sizeof(double));
 funnel = (int*) R_alloc(*nfounders, sizeof(int));
 geni = (int*) R_alloc(3, sizeof(int));
 hapi = (int*) R_alloc(3, sizeof(int));
 cf1 = (int*) R_alloc(*nfounders, sizeof(int));
 cf2 = (int*) R_alloc(*nfounders, sizeof(int));
 cf3 = (int*) R_alloc(*nfounders, sizeof(int));
 lkhd2 = (double*) R_alloc(*npos, sizeof(double));
 lkhd = (double*) R_alloc(*npos, sizeof(double));

 for (i=0; i<*nmrk; i++)
	out[i] = 0;

 for (q=0; q<*nmrk; q++)
 {
   for (i=0; i<np; i++)
	lkhd[i] = lkhd2[i] = lkhd2[i+np] = lod[q*(np)+i] = 0;

   for (i=0; i<*npos; i++)
   for (j=0; j<*nfinals; j++) 
   {
    lk = (left[i]-1)*(left[i]>0) + q*(left[i]<0);
    rk = (right[i]-1)*(right[i]>0) + q*(right[i]<0);
    mk = (mid[i]-1)*(mid[i]>0) + q*(mid[i]<0);

    // calculate haplotype probabilities, then sum over all possibilities?
    if (*nfounders==4) {
	hp4way(r12[i], r23[i], r13[i], p, x);
	for (k=0; k<10; k++)	hprob[k] = p[k]*x[k]; 
    	if ((q>0)&&(q<*npos-1)) {
		hp4way(.5, .5, r13[i], plod, xlod);
		for (k=0; k<10; k++) 	hproblod[k] = plod[k]*xlod[k];
	}
    }

    if (*nfounders==8) {
	hp8way(r12[i], r23[i], r13[i], hprob);
	if ((q>0)&&(q<*npos-1)) 
		hp8way(.5, .5, r13[i], hproblod);
    }

    geni[0] = finalg[j*(*totmrk+1)+lk+1];
    geni[1] = finalg[j*(*totmrk+1)+mk+1]; 
    geni[2] = finalg[j*(*totmrk+1)+rk+1];
    for (l=0; l<*nfounders; l++) {
 	fou[l] = founderg[l*(*totmrk)+lk];
   	fou[(*nfounders)+l] = founderg[l*(*totmrk)+mk];
	fou[2*(*nfounders)+l] = founderg[l*(*totmrk)+rk]; }

    pedfunnel(finalg[j*(*totmrk+1)], id, mother, father, funnel, *nfounders, *nped, ngen[j]);
    for (k=0; k<*nfounders; k++)
    {
 	cf1[k] = (fou[(funnel[k]-1)]==geni[0]);
	cf2[k] = (fou[*nfounders+(funnel[k]-1)]==geni[1]);
	cf3[k] = (fou[2*(*nfounders)+(funnel[k]-1)]==geni[2]);
	if (geni[0]<0) cf1[k] = 1;
	if (geni[1]<0) cf2[k] = 1;
  	if (geni[2]<0) cf3[k] = 1;
    }

   hp = 0;
   for (k=0; k<*nfounders; k++)
   for (l=0; l<*nfounders; l++)
   for (m=0; m<*nfounders; m++)
   if (cf1[k]*cf2[l]*cf3[m]) 
   {
	hapi[0] = k;
	hapi[1] = l;
	hapi[2] = m;
	hp += hprob[hclass(hapi, *nfounders)];
//	hp2 += hproblod[hclass(hapi, *nfounders)];
  }

   lkhd2[i] += log10(hp);
//   lkhd2[i+*npos] += log10(hp2);
   // end of loop over individuals
  }// end of loop over positions

//  lkhd2[*npos] = lkhd2[1+*npos];
//  lkhd2[*npos*2-1] = lkhd2[*npos*2-2];

  for (j=0; j<np; j++)
    lod[q*(np)+j] = lkhd[j] = lkhd2[j]-lkhd2[j+np];

  out[q] = whichmaxd(lkhd, np);
  }

 return;
}
