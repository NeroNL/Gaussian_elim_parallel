/*
* Unblocked LU decomposition
* based on a parallel code written by
*  Anu Rao 11/1/94 CSC 652
*
* Code modifications made by 
* Scott B. Baden, UCSD, 1/27/11
*/

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include "cblock.h"
using namespace std;

//
// Globals
//
double **A, **R;
control_block cb;

void printMat(int N,double **A){
  int i,j ;
  if (N < 10){
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            printf("%8.3g ",A[i][j]);
        }
        printf("\n"); printf("\n");
    }
    printf("\n");
 }
}

//
// External Functions
//

double getTime();
void cmdLine(int argc, char *argv[]);
void elim();
void initMat(int N, double **A, double **R, int MatCode);
void Report(int N, int NT, double t0, bool partialPivoting, int MatCode,int verif);
void ReportInput(control_block& cb);
int verify(int N);

//
// Main
//

int main(int argc,char *argv[])
{
  int i;
  double t0;

  cmdLine(argc, argv);
  printf("N = %d, nt = %d\n",cb.N,cb.NT);
  ReportInput(cb);


  /* Allocate A and  R as N by N matrices */
  assert(A = (double **)malloc(cb.N*sizeof(double *)));
  assert(R = (double **)malloc(cb.N*sizeof(double *)));
  for (i = 0; i < cb.N; i++) {
    assert(A[i] = (double *)malloc(cb.N*sizeof(double)));
    assert(R[i] = (double *)malloc(cb.N*sizeof(double)));
  } 



  initMat(cb.N,A, R,cb.MatCode);
  printMat(cb.N,A);

  //Warm up the code by running LU once, before we collect timings.
  elim();
 //
  // Don't forget to re-initialize A and R */
  initMat(cb.N, A, R, cb.MatCode);

  t0 = -getTime();
  elim();
  t0 += getTime();

  /* print results */
  printMat(cb.N,A);
  int verif = verify(cb.N);

  Report(cb.N, cb.NT,t0,cb.partialPivoting,cb.MatCode,verif);
  return(0);
}
