/*

tt_hillshade_fused.c

*/

#include "matrix.h"
#include "mex.h"
#include "topotoolbox.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
  // INPUTS: DEM, azimuth, altitude, cs
  // OUTPUTS: Hillshade
  if (nrhs != 4) {
    mexErrMsgIdAndTxt("tt3:tt_hillshade_fused:nrhs", "Four inputs required");
  }
  if (nlhs != 1) {
    mexErrMsgIdAndTxt("tt3:tt_hillshade:nlhs", "One output required");
  }

  const mxArray *demArray = prhs[0];
  float *dem = mxGetSingles(demArray);

  double azimuth = mxGetScalar(prhs[1]);
  double altitude = mxGetScalar(prhs[2]);
  double cellsize = mxGetScalar(prhs[3]);

  ptrdiff_t dims[2];
  dims[0] = mxGetM(demArray);
  dims[1] = mxGetN(demArray);

  // Null pointer checks are not needed because out-of-memory errors
  // will terminate the MEX function.
  plhs[0] = mxCreateNumericMatrix(dims[0], dims[1], mxSINGLE_CLASS, mxREAL);
  float *output = mxGetSingles(plhs[0]);

  // Converted to float to avoid compiler warnings
  hillshade_fused(output, dem, (float) azimuth, (float) altitude, (float) cellsize, dims);

}
