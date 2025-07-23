/*

tt_lowerenv.c

libtopotoolbox's lowerenv is

TOPOTOOLBOX_API
void lowerenv(float *elevation, uint8_t *knickpoints, float *distance,
              ptrdiff_t *ix, uint8_t *onenvelope, ptrdiff_t *source,
              ptrdiff_t *target, ptrdiff_t edge_count, ptrdiff_t node_count);


tt_lowerenv(z, kn, distance, source, target)
*/

#include "matrix.h"
#include "mex.h"
#include "topotoolbox.h"

#include <limits.h>
#include <stddef.h>
#include <stdint.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // Validate input and output arguments
  if (nrhs != 5) {
    mexErrMsgIdAndTxt("tt3:tt_lowerenv:nrhs", "Five inputs required");
  }
  if (nlhs != 1) {
    mexErrMsgIdAndTxt("tt3:tt_lowerenv:nlhs", "One output required");
  }

  mwSize node_count = mxGetM(prhs[0]);
  mwSize edge_count = mxGetM(prhs[3]);

  // Extract input and output array data and dimensions
  plhs[0] = mxDuplicateArray(prhs[0]);
  float *z = mxGetSingles(plhs[0]);

  uint8_t *kn = mxGetUint8s(prhs[1]);
  float *d = mxGetSingles(prhs[2]);
  ptrdiff_t *source = mxGetInt64s(prhs[3]);
  ptrdiff_t *target = mxGetInt64s(prhs[4]);

  // Allocate intermediate arrays
  ptrdiff_t *ix = mxMalloc(node_count * sizeof(*ix));
  uint8_t *onenvelope = mxMalloc(node_count * sizeof(*onenvelope));

  // Call libtopotoolbox function
  lowerenv(z, kn, d, ix, onenvelope, source, target, edge_count, node_count);

  // Destroy intermediate arrays
  mxFree(ix);
  mxFree(onenvelope);
}
