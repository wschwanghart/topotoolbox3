/*

tt_excesstopography_fmm3d.c

Computes the excess topography based on threshold slopes derived from
a three dimensional lithology using the fast marching method.
  
excess = tt_excesstopography_fmm3d(dem,lithstack,threshold_slopes,cellsize);

*/

#include "mex.h"
#include "topotoolbox.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // INPUTS: dem, lithstack, threshold_slopes, cellsize
    // OUTPUTS: excess
    if (nrhs != 4) {
        mexErrMsgIdAndTxt("ttlem3:tt_excesstopography_fmm3d:nrhs","Four inputs required");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("ttlem3:tt_excesstopography_fmm3d:nlhs","One output required");
    }

    // More validation...
    float *dem = mxGetSingles(prhs[0]);
    float *lithstack = mxGetSingles(prhs[1]);
    float *threshold_slopes = mxGetSingles(prhs[2]);
    float cellsize = (float)mxGetScalar(prhs[3]);

    mxAssert(dem,"DEM pointer is null");
    mxAssert(lithstack,"Lithstack pointer is null");
    mxAssert(threshold_slopes,"Threshold slopes pointer is null");
    
    ptrdiff_t dims[2];
    dims[0] = mxGetM(prhs[0]);
    dims[1] = mxGetN(prhs[0]);
    mwSize nlayers = mxGetM(prhs[1]);

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

    excesstopography_fmm3d(excess, heap, back, dem, lithstack, 
                           threshold_slopes, cellsize, dims, 
                           nlayers);
}
