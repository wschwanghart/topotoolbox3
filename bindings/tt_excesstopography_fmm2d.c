/*

tt_excesstopography_fmm2d.c

Computes the excess topography based on the supplied threshold slopes.
  
excess = tt_excesstopography_fmm2d(dem,threshold_slopes,cellsize);

*/

#include "mex.h"
#include "topotoolbox.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // INPUTS: dem, threshold_slopes, cellsize
    // OUTPUTS: excess
    if (nrhs != 3) {
        mexErrMsgIdAndTxt("ttlem3:tt_excesstopography_fmm2d:nrhs","Four inputs required");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("ttlem3:tt_excesstopography_fmm2d:nlhs","One output required");
    }

    // More validation...
    float *dem = mxGetSingles(prhs[0]);
    float *threshold_slopes = mxGetSingles(prhs[1]);
    float cellsize = (float)mxGetScalar(prhs[2]);

    mxAssert(dem,"DEM pointer is null");
    mxAssert(threshold_slopes,"Threshold slopes pointer is null");
    
    ptrdiff_t dims[2];
    dims[0] = mxGetM(prhs[0]);
    dims[1] = mxGetN(prhs[0]);

    mxArray *heapmat = mxCreateNumericMatrix(dims[0], dims[1], 
                                             mxINT64_CLASS, mxREAL);
    mxArray *backmat = mxCreateNumericMatrix(dims[0], dims[1], 
                                             mxINT64_CLASS, mxREAL); 

    ptrdiff_t *heap = mxGetInt64s(heapmat);
    ptrdiff_t *back = mxGetInt64s(backmat);

    mxAssert(heap,"Heap pointer is null");
    mxAssert(back,"back pointer is null");

    plhs[0] = mxCreateNumericMatrix(dims[0], dims[1], 
                                    mxSINGLE_CLASS, mxREAL);
    float *excess = mxGetSingles(plhs[0]);
    mxAssert(excess,"Excess pointer is null");

    excesstopography_fmm2d(excess, heap, back, dem, 
                           threshold_slopes, cellsize, dims);
}
