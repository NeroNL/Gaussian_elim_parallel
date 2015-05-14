/*
* Unblocked LU decomposition
* based on a parallel code written by
*  Anu Rao 11/1/94
*
* Code modifications made by 
* Scott B. Baden, UCSD, 5/3/15
*/

#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <thread>
#include "barrier.h"
#include "cblock.h"
using namespace std;

//
// Globals
//

extern double **A, **R;
extern control_block cb;
barrier count(cb.NT);
int tmp_count = 0;

//
// External Functions
//
double getTime();


//
// Gaussian Elmination
//

void serial_elim(){
  int i, j, k, Mx;

  // If we get stuck, we can count the number of row swaps
  // Do this for a small matrix
  // int swaps = 0;

  for ( k = 0; k < cb.N; k++ ) {
      if (cb.partialPivoting){ /* Partial Pivoting */
        Mx = k;
        for ( i = k+1; i < cb.N; i++ ) {
            if (fabs(A[i][k]) > fabs(A[Mx][k]))
                Mx = i;
        }
        if (Mx > k){
    //        swaps++;
            for ( j = k; j < cb.N; j++ ){
                double t = A[Mx][j];
                A[Mx][j] = A[k][j];
                A[k][j] = t;
            }
        }
      } /* End Partial Pivoting */

    for ( i = k+1; i < cb.N; i++ ) 
      A[i][k] /= A[k][k];  

    for ( i = k+1; i < cb.N; i++ ) {
      double Aik = A[i][k];
      double *Ai = A[i];
      for ( j = k+1; j < cb.N; j++ ) 
        Ai[j] -= Aik * A[k][j];
    }  
  }
}


void parallel_elim(int startIndex, int increment, int k0){
    int k = k0;
    while(k < cb.N){
      for ( int i = startIndex+k+1; i < cb.N; i+=increment ) {
        A[i][k] /= A[k][k];
      }

      count.bsync(startIndex);
      
      for ( int i = startIndex+k+1; i < cb.N; i+=increment ) {
        double Aik = A[i][k];
        double *Ai = A[i];
        for ( int j = k+1; j < cb.N; j++ ) 
          Ai[j] -= Aik * A[k][j];
      }

      count.bsync(startIndex);
      ++k;
    }
    count.bsync(startIndex);
}

void partialPivoting_parallel(int k, int Mx){
  for ( int j = k; j < cb.N; j++ ){
      double t = A[Mx][j];
      A[Mx][j] = A[k][j];
      A[k][j] = t;
  }
}

void elim(){
 if(cb.NT == 1){
    serial_elim();
 } 
 else{
      int i = 0, k = 0;
      thread* thrd = new thread[(cb.NT-1)];
  
      //cyclic partitioning
      //for(k = 0; k < cb.N; k++){
        //if(k == 0){
          for( i = 0; i < cb.NT-1; i++){
            thrd[i] = thread(parallel_elim, i, cb.NT, k);
          }
        //}
      //}
        parallel_elim(i, cb.NT, k);

        for(int i = 0; i < cb.NT-1; i++)
            thrd[i].join();

    }
}
