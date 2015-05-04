//
// Do not change the code in this file, as doing so
// could cause your submission to be graded incorrectly
//
// Process command line arguments
// 
#include <assert.h>
#include <getopt.h>
#include <stdlib.h>
#include <iostream>
#include "lu.h"
#include "cblock.h"
using namespace std;

extern control_block cb;
// Process command line arguments
void cmdLine(int argc, char *argv[]) {
 static struct option long_options[] = {
        {"n", required_argument, 0, 'n'},
        {"nt", required_argument, 0, 't'},
        {"p", no_argument, 0, 'p'},
        {"l", no_argument, 0, 'l'},
        {"h", no_argument, 0, 'h'},

 };

 // Default values
 cb.N = 16; cb.NT  = 1; cb.MatCode = D_MAT; cb.partialPivoting = false;
 int ac;
 for(ac=1;ac<argc;ac++) {
    int c;
    while ((c=getopt_long(argc,argv,"n:t:plh",long_options,NULL)) != -1){
        switch (c) {

	    // Size of the computational box
            case 'n':
                cb.N = atoi(optarg);
                break;

	    // # threads
            case 't':
                cb.NT = atoi(optarg);
                break;

	    // Use partial pivoting
            case 'p':
                cb.partialPivoting = true;
                break;

	    // Use the 'Laplacian' matrix (5-point Laplacian in 2D)
            // Doesn't need partial pivoting
            case 'l':
                cb.MatCode = L_MAT;
                break;

	    // Use a Hadamard matrix 
            // Requires partial pivoting
            case 'h':
                cb.MatCode = H_MAT;
                break;

	    // Error
            default:
                cout << "Usage: lu [-n <domain size>] [-t <# threads>] [-l <use 5-point Laplacian>[-p do partial pivoting] [-h <use Hadamard Matrix>]";
                cout << endl;
                exit(-1);
            }
    }
 }
}
