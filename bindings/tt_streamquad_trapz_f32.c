/*

tt_streamquad_trapz_f32.c

libtopotoolbox's lowerenv is

TOPOTOOLBOX_API
void streamquad_trapz_f32(float *integral, float *integrand, 
              ptrdiff_t *source, ptrdiff_t *target, float *weight,
              ptrdiff_t edge_count);


tt_streamquad_trapz_f32(integral, integrand, source, target, weight)
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
    mexErrMsgIdAndTxt("tt3:tt_streamquad_trapz32:nrhs", "Four inputs required");
  }
  if (nlhs != 1) {
    mexErrMsgIdAndTxt("tt3:tt_streamquad_trapz32:nlhs", "One output required");
  }

  mwSize edge_count = mxGetM(prhs[1]);
  mwSize r = mxGetM(prhs[0]);
  mwSize c = mxGetN(prhs[0]);

  // Extract input and output array data and dimensions
  plhs[0] = mxCreateNumericMatrix(r,c,mxSINGLE_CLASS, mxREAL);
  float *integral = mxGetSingles(plhs[0]);
  float *integrand = mxGetSingles(prhs[0]);
  ptrdiff_t *source = mxGetInt64s(prhs[1]);
  ptrdiff_t *target = mxGetInt64s(prhs[2]);
  float *w = mxGetSingles(prhs[3]);

  // Call libtopotoolbox function
  streamquad_trapz_f32(integral, integrand, source, target, w, edge_count);

}
