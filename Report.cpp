// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
#include <stdlib.h>
#include <iostream>
#include "cblock.h"

using namespace std;

/* Computes the gigaflops rate  */
double gflops(int n, double t){
    int64_t n3 = n*n; n3 *= n;
    int64_t flops = (n3*2)/3; 
    double flop_rate = (double) flops / t;
    return ( flop_rate/1.0e9);
}

void Report(int N, int NT, double tP, bool partialPivoting, int MatCode,int verif)
{
  double gf = gflops(N,  tP);
  char piv = partialPivoting ? 'Y' : 'N';
  printf("N=%d: %f seconds (%f gflops)\n\n",N,tP,gf);
  printf("    N, NT, tP, gflops -- Partial Pivoting? Verified? MtxCode\n");
  printf("@> %d %d %f %f %1c %1c %1c\n\n",N,NT,tP,gf,piv,verif,MatCode);
}
void ReportInput(control_block& cb){

  switch (cb.MatCode){
      case 'L':
          cout << "Using the '5-point' Matrix" << endl;
          break;

      case 'H':
          cout << "Using a 'Hadamard' Matrix" << endl;
          break;

      case 'D':
          cout << "Using a simple matrix 2*I" << endl;
          break;

      default:
          cerr << "No such matrix." << endl;
          break;
  }

  if (cb.partialPivoting)
      cout << "Using partial pivoting\n";
  else
      cout << "NOT Using partial pivoting\n";

}
