/*

tt_has_topotoolbox.c

Returns true

This is used to test the compilation and linking procedure for
libtopotoolbox
*/

#include "mex.h"
#include "topotoolbox.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    // INPUTS: -
    // OUTPUTS: flag
    if (nrhs != 0) {
        mexErrMsgIdAndTxt("tt3:tt_has_topotoolbox:nrhs","Zero inputs required");
    }
    if (nlhs != 1) {
        mexErrMsgIdAndTxt("tt3:tt_has_topotoolbox:nlhs","One output required");
    }

    const mwSize dims[1] = {1};
    plhs[0] = mxCreateLogicalArray(1,dims);
    mxLogical *flag = mxGetLogicals(plhs[0]);
    mxAssert(flag,"Flag pointer is null");

    *flag = has_topotoolbox();
}
