/*

tt_gwdt_computecosts.c

COSTS = tt_gwdt_computecosts(DEMF, DEM, FLATS)

*/

#include "matrix.h"
#include "mex.h"
#include "topotoolbox.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {  
  // Validate input and output arguments
  if (nrhs != 3) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt_computecosts:nrhs", "Three inputs required");
  }
  if (nlhs != 1) {
    mexErrMsgIdAndTxt("tt3:tt_gwdt_computecosts:nlhs", "One output required");
  }
  
  // Extract input and output array data and dimensions

  // Allocate necessary intermediate arrays

  // Call libtopotoolbox function
  gwdt_computecosts(costs, conncomps, flats, original_dem, filled_dem, dims);

  // Destroy intermediate arrays
}
