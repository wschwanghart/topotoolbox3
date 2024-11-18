/*

tt_gwdt_computecosts.c

*/

#include "matrix.h"
#include "mex.h"
#include "topotoolbox.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Validate input and output arguments

  // Extract input and output array data and dimensions

  // Allocate necessary intermediate arrays

  // Call libtopotoolbox function
  gwdt_computecosts(costs, conncomps, flats, original_dem, filled_dem, dims);

  // Destroy intermediate arrays
}
