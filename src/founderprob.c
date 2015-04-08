#include "mpMap.h"

// assumes unique value match
// if more than one match returns the last
int which(int x, int *vec, int length)
{
  int i;
  int whnum=-1;
  for (i=0; i<length; i++)
	if (vec[i]==x) whnum=i;
  return (whnum);
}

int sumi(int *x, int length)
{
  int i, sum=0;
 
  for (i=0; i<length; i++)
	sum += x[i];

  return(sum);
}

double sumd(double *x, int length)
{ 
  int i;
  double sum=0;

  for (i=0; i<length; i++)
	  sum += x[i];

  return(sum);
}

int whichmaxd(double *vec, int length)
{
  int i, wm=0;
  double x;

  x=vec[0];
  for (i=1; i<length; i++)
  if (vec[i]>x)
  {
    x=vec[i];
    wm = i;
  }

  return(wm);
}
// if there are ties, returns the first
int whichmaxi(int *vec, int length)
{
  int i, wm=0;
  int x;

  x=vec[0];
  for (i=1; i<length; i++)
  if (vec[i]>x)
  {
    x=vec[i];
    wm = i;
  }
  return(wm);
}



int maxi(int *vec, int length)
{
  int i;
  int x=vec[0];

  for (i=1; i<length; i++)
  if (vec[i]>x)
	x=vec[i];

  return(x);
}


double maxd(double *vec, int length)
{
  int i;
  double x=vec[0];

  for (i=1; i<length; i++)
  if (vec[i]>x)
    x=vec[i];

  return(x);
}

// Note: 4 way probabilities are actually x*p
void hp4way(double r12, double r23, double r13, double *p, double *x)
{
 int i;

 p[0] = .25*(1-r12+1-r13-r23);
 p[1] = .25*(1-r12+r13-1+r23);
 p[2] = .25*(r12-r23+r13);
 p[3] = .25*(r12-r13+r23);
 p[4] = .25*(1-r12);
 p[5] = .25*(1-r23);
 p[6] = r12/4;
 p[7] = r23/4;
 p[8] = r13/4;
 p[9] = (1-r13)/4;

 for (i=0; i<4; i++)
   x[i] = .25*(1/(1+2*r12)+1/(1+2*r13)-2*r23/(1+2*r23));

 x[4] = x[6] = .25*(1/(1+2*r12)+2*r13/(1+2*r13)-1/(1+2*r23));
 x[5] = x[7] = .25*(2*r12/(1+2*r12)+1/(1+2*r23)-1/(1+2*r13));
 x[8] = x[9] = .25*(2*r12/(1+2*r12)+1/(1+2*r13)-1/(1+2*r23));

 return;
}

void hp8way(double r12, double r23, double r13, double *hp)
{
  double *p, *x;

  x = (double*)R_alloc(10, sizeof(double));
  p = (double*)R_alloc(10, sizeof(double));

  hp4way(r12, r23, r13, p, x);

  hp[0] = x[0]*p[0]*p[0];
  hp[1] = x[0]*p[0]*p[1];
  hp[2] = x[0]*p[0]*p[2];
  hp[3] = x[1]*p[1]*p[4];
  hp[4] = x[2]*p[2]*p[5];
  hp[5] = x[4]*p[4]*p[4];
  hp[6] = x[5]*p[5]*p[5];
  hp[7] = x[0]*p[0]*p[3];
  hp[8] = x[3]*p[3]*p[9];
  hp[9] = x[9]*p[9]*p[9];
  hp[10] = x[1]*p[1]*p[6];
  hp[11] = x[4]*p[4]*p[6];
  hp[12] = x[2]*p[2]*p[7];
  hp[13] = x[6]*p[6]/8;
  hp[14] = x[5]*p[5]*p[7];
  hp[15] = x[8]*p[8]/8;
  hp[16] = x[7]*p[7]/8;
  hp[17] = x[3]*p[3]*p[8];
  hp[18] = x[9]*p[9]*p[8];

  return;
}

void ftobin(int *hap, int x, int n)
{
  int i;
  int r=x;

  hap[0] = (r >= (1 << (n-1)));
  for (i=1; i<n; i++)
  {
    r -= hap[i-1]*(1 << (n-i));
    hap[i] = (r >= (1 << (n-i-1)));
  }

  return;
}




int hclass(int *hap, int nfounders)
{
  int n = log(nfounders)/log(2);
  int i, j;
  int h1bin[n], h2bin[n], neq[n];
  int vec[3];
  int pair1[3]={0, 0, 1}, pair2[3]={1, 2, 2};
  int score;
  int scorevec[n];
  int match4[10] = {0,4,12,10,8,24,17,25,23,20};
  int match8[19] = {0,5,20,10,40,15,60,17,34,51,26,31,41,47,61,59,62,38,55};
  int hc=0;

  for (i=0; i<3; i++)
  {
    ftobin(h1bin, hap[pair1[i]], n);
    ftobin(h2bin, hap[pair2[i]], n);
    for (j=0; j<n; j++)
	    neq[j] = (h1bin[j]!=h2bin[j]);
    vec[i] = (n-whichmaxi(neq, n))*(maxi(neq, n)>0);
    scorevec[i] = pow(n+1, 2-i)*vec[i];
  }

  score = sumi(scorevec, 3);

  if (nfounders==4) hc = which(score, match4, 10);
  if (nfounders==8) hc = which(score, match8, 19);

  return(hc);
}

void pedfunnel(int i, int *id, int *mother, int *father, int *ofunnel, int nfounders, int nped, int ngen)
{
  int subj = which(i, id, nped);
  
  int m, f, mm, mf, ff, fm, mmm, mmf, mfm, mff, ffm, fff, fmm, fmf;

  // put in error check if which function returns a negative number

  // to remove selfing
  while(mother[subj]==father[subj])
	subj = which(mother[subj], id, nped);
  int aiCounter = ngen;
  while(aiCounter)
  {
	subj = which(mother[subj], id, nped);
	aiCounter--;
  }
  m = mother[subj];
  f = father[subj];

  mf = father[which(m, id, nped)];
  mm = mother[which(m, id, nped)];

  ff = father[which(f, id, nped)];
  fm = mother[which(f, id, nped)];


 if (nfounders==8) {
  ofunnel[0] = mmm = mother[which(mm, id, nped)];
  ofunnel[1] = mmf = father[which(mm, id, nped)];
  ofunnel[2] = mfm = mother[which(mf, id, nped)];
  ofunnel[3] = mff = father[which(mf, id, nped)];

  ofunnel[6] = ffm = mother[which(ff, id, nped)];
  ofunnel[7] = fff = father[which(ff, id, nped)];
  ofunnel[4] = fmm = mother[which(fm, id, nped)];
  ofunnel[5] = fmf = father[which(fm, id, nped)];
  }

 if (nfounders==4)
 {
	ofunnel[0] = mm;
	ofunnel[1] = mf;
	ofunnel[2] = fm;
	ofunnel[3] = ff;
 }
	
 return;
}

// Note that this requires that the id, mother, father portions of the pedigree are a) integers and b) the founders are in the first 4/8 positions
// what else needs to be input? step size? some way of scanning the interval
// also need framework map..
//
// mrkstar is the marker to be positioned; mrk1 and mrk 2 are the additional markers at each point determined by r12 and r23 (assume no interference)
void hap3ptfull(int *finalg, int *founderg, int *id, int *mother, int *father, int *nfinals, int *nfounders, int *nmrk, int *nped, int* ngen, double *out, double *r12, double *r23, double *r13, int *left, int *mid, int *right, int *npos)
{
 int i, j, k, l, m;
 int *funnel, *geni, *hapi;
 double *hprob, *p, *x;
 int *cf1, *cf2, *cf3;
 int *fou;
 double *hp;
 double sum;

 hp = (double*) R_alloc(*nfounders, sizeof(double));
 fou = (int*) R_alloc(*nfounders*3, sizeof(int));
 p = (double*) R_alloc(10, sizeof(double));
 x = (double*) R_alloc(10, sizeof(double));
 hprob = (double*) R_alloc(20, sizeof(double));
 funnel = (int*) R_alloc(*nfounders, sizeof(int));
 geni = (int*) R_alloc(*nmrk, sizeof(int));
 hapi = (int*) R_alloc(*nmrk, sizeof(int));
 cf1 = (int*) R_alloc(*nfounders, sizeof(int));
 cf2 = (int*) R_alloc(*nfounders, sizeof(int));
 cf3 = (int*) R_alloc(*nfounders, sizeof(int));
 

 // loop over positions
 for (i=0; i<*npos; i++)
 {
   // calculate haplotype probabilities, then sum over all possibilities?
   if (*nfounders==4) {
	hp4way(r12[i], r23[i], r13[i], p, x);
	for (j=0; j<10; j++)	hprob[j] = p[j]*x[j]; }

   if (*nfounders==8)
	hp8way(r12[i], r23[i], r13[i], hprob);

  for (j=0; j<*nfinals; j++)
   {
     for (k=0; k<*nfounders; k++)
     out[j*(*npos)*(*nfounders)+i*(*nfounders)+k] = 0;
 
     geni[0] = finalg[j*(*nmrk+1)+left[i]];
     geni[1] = -(mid[i]<=0)+(mid[i]>0)*(finalg[j*(*nmrk+1)+mid[i]]); 
     geni[2] = finalg[j*(*nmrk+1)+right[i]];
     for (l=0; l<*nfounders; l++) {
	fou[l] = founderg[l*(*nmrk)+left[i]-1];
   	fou[(*nfounders)+l] = (mid[i]>0)*(founderg[l*(*nmrk)+mid[i]-1]);
	fou[2*(*nfounders)+l] = founderg[l*(*nmrk)+right[i]-1]; }

     pedfunnel(finalg[j*(*nmrk+1)], id, mother, father, funnel, *nfounders, *nped, ngen[j]);
 
     for (k=0; k<*nfounders; k++)
     {
 	cf1[k] = (fou[(funnel[k]-1)]==geni[0]);
	cf2[k] = (fou[*nfounders+(funnel[k]-1)]==geni[1]);
	cf3[k] = (fou[2*(*nfounders)+(funnel[k]-1)]==geni[2]);
	if (geni[0]<0) cf1[k] = 1;
	if (geni[1]<0) cf2[k] = 1;
 	if (geni[2]<0) cf3[k] = 1;
     }
	
    for (k=0; k<*nfounders; k++)
     hp[k] = 0;

    for (k=0; k<*nfounders; k++)
    for (l=0; l<*nfounders; l++)
    for (m=0; m<*nfounders; m++)
    if (cf1[k]*cf2[l]*cf3[m]) 
    {
	hapi[0] = k;
	hapi[1] = l;
	hapi[2] = m;
	// check hclass function with inputs... how does ti change with 3 loci
	if (i==0)
		hp[k] += hprob[hclass(hapi, *nfounders)];
	else if (i==(*npos-1))
		hp[m] += hprob[hclass(hapi, *nfounders)];
	else
		hp[l] += hprob[hclass(hapi, *nfounders)];
    }

	sum = 0; 
	for(k = 0; k < *nfounders; k++)
	{
		sum += hp[k];
	}

    for (k=0; k<*nfounders; k++)
    out[j*(*npos)*(*nfounders)+i*(*nfounders)+k] = hp[k]/sum;

   } // end of loop over individuals
 }// end of loop over positions

 return;
}
