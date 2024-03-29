//
// include necessary system headers
//
#include <cmath>
#include <mex.h>
#include <array>

#include <iostream>

/* MEX gateway */

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray *prhs[] )
     
{ 
    // In: diagonal (nc x m), N (nf x 2)
    // Out: Face diagonal of (2*nf x m)
    if (nrhs != 2) { 
	    mexErrMsgTxt("2 input arguments required."); 
    } else if (nlhs > 1) {
	    mexErrMsgTxt("Wrong number of output arguments."); 
    } 
    double * diagonal = mxGetPr(prhs[0]);
    double * N = mxGetPr(prhs[1]);

    int nc = mxGetM(prhs[0]);
    int m = mxGetN(prhs[0]);
    int nf = mxGetM(prhs[1]);
        
    plhs[0] = mxCreateDoubleMatrix(2*nf, m, mxREAL);
    double * result = mxGetPr(plhs[0]);
    #ifdef _WIN32
        #pragma omp parallel for
    #else
        #pragma omp parallel for collapse(2)
    #endif
    for(int j=0;j<m;j++){
        for(int i=0;i<2*nf;i++){
            int cell_inx = N[i]-1;
            result[j*2*nf + i] = diagonal[nc*j + cell_inx]/2.0;
        }
    }
    return;
}


