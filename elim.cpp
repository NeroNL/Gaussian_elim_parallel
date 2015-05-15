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


void parallel_elim(int startIndex, int increment){
    int k = 0, Mx = 0;
    int i = 0, j = 0;

    while(k < cb.N){
      /* Partial Pivoting */
      /*if (cb.partialPivoting){ 
        Mx = k;
        for ( int i = startIndex+k+1; i < cb.N; i+=increment ) {
            if (fabs(A[i][k]) > fabs(A[Mx][k]))
                Mx = i;
        }
        partial_index[startIndex] = Mx;
        count.bsync(startIndex);
        partialPivoting_parallel(k, startIndex, increment);
        count.bsync(startIndex);
      }*/ /* End Partial Pivoting */

      for ( i = startIndex+k+1; i < cb.N; i+=increment ) {
        A[i][k] /= A[k][k];
        if(i == 3)
          cout << "stupid mistake " << startIndex << endl;
      }

      //count.bsync(startIndex);

      for ( i = startIndex+k+1; i < cb.N; i+=increment ) {
        double Aik = A[i][k];
        double *Ai = A[i];
        for ( j = k+1; j < cb.N; j++ ) 
          Ai[j] -= Aik * A[k][j];
      }

      ++k;
      count.bsync(startIndex);

    }

    //cout << "k is " << k << " start index is " << startIndex << endl;
    count.bsync(startIndex);
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
            thrd[i] = thread(parallel_elim, i, cb.NT);
          }
        //}
      //}
        parallel_elim(i, cb.NT);

        for(int i = 0; i < cb.NT-1; i++)
            thrd[i].join();

    }
}
