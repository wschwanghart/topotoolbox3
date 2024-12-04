/*

tt_gwdt_computecosts.c

COSTS = tt_gwdt_computecosts(DEM, DEMF, FLATS)

*/

#include "matrix.h"
#include "mex.h"
#include "topotoolbox.h"

#include <limits.h>
#include <stddef.h>
#include <stdint.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Validate input and output arguments
  if (nrhs != 3) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt_computecosts:nrhs", "Three inputs required");
  }
  if (nlhs != 1) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt_computecosts:nlhs", "One output required");
  }
  // Extract input and output array data and dimensions
  const mxArray *demArray = prhs[0]; // DEM
  float *original_dem = mxGetSingles(demArray);

  const mxArray *demfArray = prhs[1]; // DEMF
  float *filled_dem = mxGetSingles(demfArray);

  const mxArray *flatsArray = prhs[2]; // FLATS
  int32_t *flats = mxGetInt32s(flatsArray);

  // Create the dimensions.
  //
  // This is always an array of two ptrdiff_t.
  //
  // A ptrdiff_t is equivalent to mwSignedIndex. MATLAB's mxGetM and
  // mxGetN return mwSize, which is a size_t, i.e. an unsigned
  // integer. Converting from size_t to a ptrdiff_t is unsafe, because
  // the latter is not guaranteed to be large enough to hold a
  // size_t. In general we don't check this because it is not feasible
  // to have an array of single-precision floating point numbers with
  // a size that fits in a size_t but not a ptrdiff_t. Here, I do show
  // how to check it for completeness.

  ptrdiff_t dims[2];
  mwSize m = mxGetM(demArray);
  mwSize n = mxGetN(demArray);

  if (m > PTRDIFF_MAX) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt_computecosts:mxGetM",
                      "Dimension does not fit in a ptrdiff_t");
  }
  if (n > PTRDIFF_MAX) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt_computecosts:mxGetN",
                      "Dimension does not fit in a ptrdiff_t");
  }

  dims[0] = m;
  dims[1] = n;

  // Create output array
  plhs[0] = mxCreateNumericMatrix(dims[0], dims[1], mxSINGLE_CLASS, mxREAL);
  float *costs = mxGetSingles(plhs[0]);

  // Allocate necessary intermediate arrays
  if (dims[0] > 0 && dims[1] > PTRDIFF_MAX / dims[0]) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt_computecosts:mxMalloc",
                      "Element count overflows");
  }
  ptrdiff_t count = dims[0] * dims[1];

  if (count > PTRDIFF_MAX / sizeof(ptrdiff_t)) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt_computecosts:mxMalloc",
                      "Intermediate array size overflows");
  }
  ptrdiff_t *conncomps = mxMalloc(sizeof(ptrdiff_t) * count);

  // Call libtopotoolbox function
  gwdt_computecosts(costs, conncomps, flats, original_dem, filled_dem, dims);

  // Destroy intermediate arrays
  mxFree(conncomps);
}
