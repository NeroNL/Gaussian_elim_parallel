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
int *partial_index = new int[cb.NT];

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
    if(k == 1)
      break;
  }
}



void partialPivoting_parallel(int k, int startIndex, int increment){
  int tmp = partial_index[0];
  for(int i = 0; i < cb.NT;i++){
    if(tmp < partial_index[i]){
      tmp = partial_index[i];
    }
  }

  if(tmp > k){
    for(int j = startIndex+k; j < cb.N; j+=increment){
      double t = A[tmp][j];
      A[tmp][j] = A[k][j];
      A[k][j] = t;
    }
  }
}


void parallel_elim(int TID, int NumThread){
    int k = 0, Mx = 0;
    int i = 0, j = 0;

    while(k < cb.N){
      int distance = cb.N-k-1;
      int start = (TID * distance / NumThread)+k;
      int end = ((TID+1) * distance / NumThread)+k;

      if(start < cb.N && end < cb.N){
        for ( i = start+1; i <= end; i++ ) {
          A[i][k] /= A[k][k];
          //if(k == 1)
            cout << "i is " << i << " TID is " << TID << " k is " << k << endl;
        }
      }
      i = 0;
      count.bsync(start);

      for ( i = start+1; i <= end; i++ ) {
        double Aik = A[i][k];
        double *Ai = A[i];
        if(k == 1 && i == 2 && TID == 0)
            cout << "stupid" << endl;
        for ( j = k+1; j < cb.N; j++ ) 
          Ai[j] -= Aik * A[k][j];
      }
  
    if(k == 1)
      break;
      ++k;
      count.bsync(start);

    }

    //cout << "k is " << k << " start index is " << startIndex << endl;
    //count.bsync(start);
}



void elim(){
 if(cb.NT == 1){
    serial_elim();
 } 
 else{
      int i = 0, k = 0;
      thread* thrd = new thread[(cb.NT)];
  
      //cyclic partitioning
      //for(k = 0; k < cb.N; k++){
        //if(k == 0){
          for( i = 0; i < cb.NT; i++){
            //int start = i*(cb.NT/cb.N);
            //int end = (i+1)*(cb.NT/cb.N);
            thrd[i] = thread(parallel_elim, i, cb.NT);
          }
        //}
      //}

        for(int i = 0; i < cb.NT; i++)
            thrd[i].join();

    }
}
