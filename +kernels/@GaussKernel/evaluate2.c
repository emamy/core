#include <mex.h>
#include <math.h>
/*#include <omp.h>*/

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
   double *x, *y, *res, gamma, hlp;

   int cur, i, j, k, l, m, n, n2;

   /* const mxArray *net; */

   /*
   if (nrhs != 3)
   {
      mexErrMsgTxt("Wrong number input arguments.");
   }
   else if (nlhs > 1)
   {
      mexErrMsgTxt("Too many output arguments.");
   }*/

   /*
   if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
   {
      mexErrMsgTxt("x1 must be a double matrix.");
   }
   */
   n  = mxGetM(prhs[1]);
   m  = mxGetN(prhs[1]);
   x = mxGetPr(prhs[1]);

   /*
   if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
   {
      mexErrMsgTxt("x2 must be a double matrix.");
   }
   */
   int sel = 2;
   if (nrhs != 3) sel = 1;
   n2  = mxGetM(prhs[sel]);
   k  = mxGetN(prhs[sel]);
   y = mxGetPr(prhs[sel]);
   
   if (n != n2)
    mexErrMsgTxt("GaussKernel.evaluate (mex): Array dimensions must match");
   
   /* Output */
   plhs[0] = mxCreateDoubleMatrix(m, k, mxREAL);
   res = mxGetPr(plhs[0]);

   
   /*char temp[100];*/
   /* Matrix computation */
   gamma = mxGetScalar(mxGetProperty(prhs[0],0,"Gamma"));
   
   /*#pragma omp parallel private(i,j,l,cur,hlp) shared(res,gamma,x,y,m,k)*/
   for (i = 0; i < m; i++)
   {
      for (j = 0; j < k; j++)
      {
         cur = i+j*m;
         /*sprintf(temp,"%.12f",cur);
         mexPrintf(temp);*/
         for (l = 0; l < n; l++)
         {
            hlp = x[i*n+l] - y[j*n+l];
            res[cur] += hlp*hlp;
         }
         res[cur] = exp(-gamma*res[cur]);
      }
   }
   
   /* Exponate */
   
   /*for (i = 0; i < m*n; i++)
   {
      
   }*/
}
