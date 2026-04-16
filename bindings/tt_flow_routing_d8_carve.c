/*

  [SOURCE, TARGET, COUNT] = tt_flow_routing_d8_carve(DEMF, DIST, FLATS);
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
    mexErrMsgIdAndTxt("tt3:tt_flow_routing_d8_carve:nrhs", "Three inputs required");
  }
  if (nlhs != 3) {
    mexErrMsgIdAndTxt("tt3:tt_flow_routing_d8_carve:nlhs", "Three outputs required");
  }

  float *demf = mxGetSingles(prhs[0]);
  float *costs = mxGetSingles(prhs[1]);
  int32_t *flats = mxGetInt32s(prhs[2]);

  ptrdiff_t dims[2] = {mxGetM(prhs[0]), mxGetN(prhs[0])};

  plhs[0] = mxCreateNumericMatrix(dims[0] * dims[1], 1, mxINT64_CLASS, mxREAL);
  int64_t *source = mxGetInt64s(plhs[0]);

  plhs[1] = mxCreateNumericMatrix(dims[0] * dims[1], 1, mxINT64_CLASS, mxREAL);
  int64_t *target = mxGetInt64s(plhs[1]);
  
  ptrdiff_t *node = mxMalloc(sizeof(ptrdiff_t) * dims[0] * dims[1]);
  uint8_t *direction = mxMalloc(sizeof(uint8_t) * dims[0] * dims[1]);

  flow_routing_d8_carve(node, direction, demf, costs, flats, dims, 0);

  ptrdiff_t edge_count =
      flow_routing_d8_edgelist(source, target, node, direction, dims, 0);

  plhs[2] = mxCreateDoubleScalar(edge_count);

  mxFree(node);
  mxFree(direction);
}
