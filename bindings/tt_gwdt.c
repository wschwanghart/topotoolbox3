/*

tt_gwdt.c

[AUXTOPO] = tt_gwdt(DEM, DEMF, FLATS, COSTS)

FLATS has to be the 32 bit encoded array returned by tt_identifyflats.

*/

#include "matrix.h"
#include "mex.h"
#include "topotoolbox.h"

#include <limits.h>
#include <stddef.h>
#include <stdint.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Validate input and output arguments
  if (nrhs != 4) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt:nrhs", "Four inputs required");
  }
  if (nlhs != 1) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt:nlhs", "One outputs required");
  }
  // Extract input and output array data and dimensions
  float *dem = mxGetSingles(prhs[0]);
  float *demf = mxGetSingles(prhs[1]);
  int32_t *flats = mxGetInt32s(prhs[2]);
  float *costs = mxGetSingles(prhs[3]);

  // Create the dimensions.
  ptrdiff_t dims[2] = {mxGetM(prhs[0]), mxGetN(prhs[0])};

  // Create output and intermediate arrays
  plhs[0] = mxCreateNumericMatrix(dims[0], dims[1], mxSINGLE_CLASS, mxREAL);
  float *auxtopo = mxGetSingles(plhs[0]);

  // Allocate necessary intermediate arrays
  if (dims[0] > 0 && dims[1] > PTRDIFF_MAX / dims[0]) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt:mxMalloc",
                      "Element count overflows");
  }
  ptrdiff_t count = dims[0] * dims[1];

  if (count > PTRDIFF_MAX / sizeof(ptrdiff_t)) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt:mxMalloc",
                      "Intermediate array size overflows");
  }
  
  ptrdiff_t *heap = mxMalloc(sizeof(ptrdiff_t) * count);
  ptrdiff_t *back = mxMalloc(sizeof(ptrdiff_t) * count);

  gwdt(auxtopo, NULL, costs, flats, heap, back, dims);

  // Destroy intermediate arrays
  mxFree(heap);
  mxFree(back);
}
