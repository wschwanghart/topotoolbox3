/*

tt_graphflood.c

*/

#include "matrix.h"
#include "mex.h"
#include "topotoolbox.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // INPUTS: dem, hw, BCs, Precipitations, manning, dt, dx, SFD, D8, N_iterations, step
  // OUTPUTS: hw
  if (nrhs != 11) {
    mexErrMsgIdAndTxt("tt3:tt_graphflood:nrhs", "Eleven inputs required");
  }
  if (nlhs != 1) {
    mexErrMsgIdAndTxt("tt3:tt_graphflood:nlhs", "One output required");
  }

  // DEM array
  float *dem = mxGetSingles(prhs[0]);

  // Duplicate water height array so that it serves as output
  plhs[0] = mxDuplicateArray(prhs[1]);
  float *hwout = mxGetSingles(plhs[0]);

  // Boundary conditions
  uint8_t *BCs = mxGetUint8s(prhs[2]);
  
  // Precip
  float *Precipitations = mxGetSingles(prhs[3]);

  // Manning
  float *manning = mxGetSingles(prhs[4]);

  // Dims
  size_t dims[2];
  dims[0] = mxGetN(prhs[0]);
  dims[1] = mxGetM(prhs[0]);

  // Time step length  
  float dt = (float) mxGetScalar(prhs[5]);
  
  // Spatial resolution
  float dx = (float) mxGetScalar(prhs[6]);

  // SFD
  bool SFD = (bool) mxGetScalar(prhs[7]);
  
  // D8
  bool D8 = (bool) mxGetScalar(prhs[8]);

  // N_iterations
  size_t N_iterations = (size_t) mxGetScalar(prhs[9]);
  
  // step
  float step = (float) mxGetScalar(prhs[10]);

  // Run graphflood
  graphflood_full(dem,  hwout, BCs, Precipitations, manning, dims, 
                  dt, dx, SFD, D8, N_iterations, step);

}
