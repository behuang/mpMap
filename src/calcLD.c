#include "mpMap.h"
const int mask[8][8] =
{
	{0, 1, 2, 2, 2, 2, 2, 2},
	{1, 0, 2, 2, 2, 2, 2, 2},
	{2, 2, 0, 1, 2, 2, 2, 2},
	{2, 2, 1, 0, 2, 2, 2, 2},
	{2, 2, 2, 2, 0, 1, 2, 2},
	{2, 2, 2, 2, 1, 0, 2, 2},
	{2, 2, 2, 2, 2, 2, 0, 1},
	{2, 2, 2, 2, 2, 2, 1, 0}
};
void pr2pt4way(double r, double *prob)
{
	prob[0] = (1-r)/(4+8*r);
	prob[1] = r/(4+8*r);
	prob[2] = r/(4+8*r);
}
void pr2pt8way(double r, double *prob)
{
	prob[0] = (1-r)*(1-r)/(8+16*r);
	prob[1] = r*(1-r)/(8+16*r);
	prob[2] = r/2/(8+16*r);
}
void pr2ptirip(double r, int s, int n, double *prob)
{
	int x = 1 << n;
	prob[0]=(pow(1-r, n+s-1)/x+(2*r+1-pow(1-r, s-1))/x/x)/(1+2*r); 
	prob[1]=prob[2]=(1-x*prob[0])/(x*(x-1));
}
void calcLD(int *finalg, int *founderg, int *id, int *mother, int *father, int *pair1, int *pair2, int *nfinals, int *nfounders, int *nmrk, int *npairs, int *nped, int *ngen, double *rpair, double *ldw, double *ldlew, double *lddelta, double *ldr2)
{
  int i, j, k, l, m;
  int p1, p2;
  int *cf1, *cf2, *funnel, genp1, genp2;
  double *probclass, delta, dmax=0, theta;
  double *tmptable, *table, *p, *q;
  int n = (int) log(*nfounders)/log(2);
  double totsum=0, wt=0;

  probclass = (double*) R_alloc(3, sizeof(double));
  table = (double*) R_alloc((*nfounders)*(*nfounders), sizeof(double));
  tmptable = (double*) R_alloc((*nfounders)*(*nfounders), sizeof(double));
  p = (double*) R_alloc(*nfounders, sizeof(double));
  q = (double*) R_alloc(*nfounders, sizeof(double));

  cf1 = (int*) R_alloc(*nfounders, sizeof(int));
  cf2 = (int*) R_alloc(*nfounders, sizeof(int));
  funnel = (int*) R_alloc(*nfounders, sizeof(int));

  for (i=0; i<(*npairs); i++)
 	ldw[i] = ldlew[i] = lddelta[i] = ldr2[i] = 0;

  for (k=0; k<*npairs; k++)
  {
    p1 = pair1[k]-1;
    p2 = pair2[k]-1;
    theta = rpair[k];

    for (i=0; i<(*nfounders)*(*nfounders); i++)
	table[i]=0;

    for (i=0; i<*nfounders; i++)
	p[i] = q[i] = 0;

    totsum=0;
    for (i=0; i<*nfinals; i++)
    {
        if (ngen[i]==0) {
           if (*nfounders==4)  pr2pt4way(theta, probclass);
           if (*nfounders==8)  pr2pt8way(theta, probclass);
        }
        if (ngen[i] > 0)
           pr2ptirip(theta, ngen[i], n, probclass);
	genp1 = finalg[i*(*nmrk+1)+p1+1];
	genp2 = finalg[i*(*nmrk+1)+p2+1];
	wt=0;

   	pedfunnel(finalg[i*(*nmrk+1)], id, mother, father, funnel, *nfounders, *nped, ngen[i]);
	for (j=0; j<(*nfounders)*(*nfounders); j++)
		tmptable[j] = 0;

    	if ((genp1>-1)*(genp2>-1))
     	{
      	  for (j=0; j<*nfounders; j++)
      	  {
	    cf1[j] = (founderg[(funnel[j]-1)*(*nmrk)+p1]==genp1);
	    cf2[j] = (founderg[(funnel[j]-1)*(*nmrk)+p2]==genp2);
      	  }

	  for (j=0; j<*nfounders; j++)
	  for (l=0; l<*nfounders; l++)
	  if (cf1[j]*cf2[l]) {
		tmptable[j*(*nfounders)+l] = probclass[mask[j][l]];
		wt += tmptable[j*(*nfounders)+l];
	  }
	 
	  for (j=0; j<(*nfounders)*(*nfounders); j++)
		table[j] += tmptable[j]/wt;
        } // end of NA check
     } // end of nfinals loop
  
   totsum = sumd(table, *nfounders*(*nfounders));

   for (i=0; i<*nfounders; i++)
   for (j=0; j<*nfounders; j++)
   {
	p[i] += table[i*(*nfounders)+j]/totsum;
	q[j] += table[i*(*nfounders)+j]/totsum;
   }

   for (i=0; i<*nfounders; i++)
   for (j=0; j<*nfounders; j++)
   {
     m = i*(*nfounders)+j;
     delta = table[m]/totsum-p[i]*q[j];
     if (delta>0)
	     dmax = q[j]*(1-p[i])*(q[j]*(1-p[i]) < (1-q[j])*p[i]) + 
		     p[i]*(1-q[j])*(p[i]*(1-q[j]) <= q[j]*(1-p[i]));
     if (delta<=0)
	     dmax = (1-p[i])*(1-q[j])*((1-p[i])*(1-q[j]) < p[i]*q[j]) +
		     p[i]*q[j]*(p[i]*q[j] <= (1-p[i])*(1-q[j]));

     ldw[k] += delta*delta/p[i]/q[j];
     ldlew[k] += p[i]*q[j]*fabs(delta)/dmax;
     lddelta[k] += delta*delta; 
	ldr2[k] += delta*delta/(1-p[i])/(1-q[j])/p[i]/q[j]*table[m]/totsum;
   }

   ldw[k] = ldw[k]/(*nfounders-1);
  } // end of npairs loop

 return;
}
