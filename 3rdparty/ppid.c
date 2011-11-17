#include <sys/types.h>
#include <unistd.h>
#include "mex.h"
/*
 * Outputs the current process ID of Matlab.
 *
 * @author Massimiliano Salsi (http://www.mathworks.com/matlabcentral/newsreader/author/100603)
 *
 * Code from http://www.mathworks.com/matlabcentral/newsreader/view_thread/164015
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray*prhs[] )
{
    plhs[0] = mxCreateDoubleScalar( getppid() );
    return;
}