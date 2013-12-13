#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include "rngs.h"

int which(int x, int *vec, int length);
int sumi(int *x, int length);
double sumd(double *x, int length);
int whichmaxd(double *vec, int length);
int whichmaxi(int *vec, int length);
int maxi(int *vec, int length);
double maxd(double *vec, int length);

void ftobin(int *hap, int x, int n);
int hclass(int *hap, int nfounders);
double kmf(double x);
double hmf(double x);
void quickSort(int *arr, int elements);

void hp8way(double r12, double r23, double r13, double *hp);
void hp4way(double r12, double r23, double r13, double *p, double *x);

void pr2pt4way(double r, double *prob);
void pr2pt8way(double r, double *prob);
void pr2ptirip(double r, int s, int n, double *prob);

extern const int mask[8][8];
void pedfunnel(int i, int *id, int *mother, int *father, int *ofunnel, int nfounders, int nped, int ngen);

void calcLD(int *finalg, int *founderg, int *id, int *mother, int *father, int *pair1, int *pair2, int *nfinals, int *nfounders, int *nmrk, int *npairs, int *nped, int *ngen, double *rpair, double *ldw, double *ldlew, double *lddelta, double *ldr2);

void hap3ptfull(int *finalg, int *founderg, int *id, int *mother, int *father, int *nfinals, int *nfounders, int *nmrk, int *nped, int* ngen, double *out, double *r12, double *r23, double *r13, int *left, int *mid, int *right, int *npos);

void pf(int *obs, int*id, int *mother, int *father, int *nfinals, int *nfounders, int *nped, int* ngen, int *out);
void creategam(double *genome, int nmrk, int *ibd, int *gamete);
void gengeno(double *genome, double *genome2, int *id, int *mother, 
		int *father, int *nmrk, int *nfdr, int *nped, int *seed, 
		int *transpos, int *transval, int *out);

