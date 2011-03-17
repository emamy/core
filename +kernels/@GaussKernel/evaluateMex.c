#include <mex.h>
#include <math.h>
#include <omp.h>

/* Implements the sum(x.^2,1) matlab command for a matrix x */
void sumsq(double* xsq, double* x, int n, int m) {
    int i, j, *v, vfirst;

#ifdef DEBUG
    char temp[100];
#endif

    #pragma omp parallel for shared(xsq,x,n,m) private(i,j)
    for (j=0;j<m;j++) {
        xsq[j] = 0;
/*
        vfirst = x[j*n]; v = &vfirst;
*/
        for (i=0;i<n;i++) {
            xsq[j] += x[j*n+i]*x[j*n+i];
/*
            xsq[j] += v[i]*v[i];
*/
#ifdef DEBUG
            sprintf(temp, "j(m)=%d(%d), i(n)=%d(%d), x=%.12f\n", j, m, i, n, x[j*n+i]);
            mexPrintf(temp);
#endif
        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    double *x, *y, *res, gamma, hlp;
    
    int i, j, k, l, m, n;
    
#ifdef DEBUG
    char temp[100];
#endif

    /* const mxArray *net; */

    /*
     * if (nrhs != 3)
     * {
     * mexErrMsgTxt("Wrong number input arguments.");
     * }
     * else if (nlhs > 1)
     * {
     * mexErrMsgTxt("Too many output arguments.");
     * }*/

    /*
     * if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]))
     * {
     * mexErrMsgTxt("x1 must be a double matrix.");
     * }
     */
    n  = mxGetM(prhs[1]);
    m  = mxGetN(prhs[1]);
    x = mxGetPr(prhs[1]);

    /*
     * if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2]))
     * {
     * mexErrMsgTxt("x2 must be a double matrix.");
     * }
     */

    double xsq[m];
    sumsq(xsq, x, n, m);

    /* Get gamma value */
    gamma = mxGetScalar(mxGetProperty(prhs[0], 0, "Gamma"));
    
#ifdef DEBUG
    sprintf(temp, "Gamma=%.12f\n", gamma);
    mexPrintf(temp);
#endif

    /* One argument case (symmetric matrix) */
    if (nrhs == 2) {
        /* Output init */
        plhs[0] = mxCreateDoubleMatrix(m, m, mxREAL);
        res = mxGetPr(plhs[0]);
        
#ifdef DEBUG
        sprintf(temp, "Call with one parameter!\n");
        mexPrintf(temp);
#endif
        /* Diagonal with ones */
        for(i=0;i<m;i++) {
            res[i*m+i] = 1;
        }

        /* Only compute lower left triangular & copy values */
        #pragma omp parallel for shared(res,xsq,x,n,m) private(hlp,i,j,l)
        for (i = 1; i < m; i++) {
            for (j = 0; j < i; j++) {
                hlp=0;
                for (l = 0; l < n; l++) {
#ifdef DEBUG
                    sprintf(temp, "ScaP i=%d, j=%d, l=%d, x=%.12f, y=%.12f\n", i, j, l, x[i*n+l], y[j*n+l]);
                    mexPrintf(temp);
#endif
                    hlp += x[i*n+l]*x[j*n+l];
                }
#ifdef DEBUG
                sprintf(temp, "xsq[i]=%.12f, xsq[j]=%.12f, x*y=%.12f\n", xsq[i], xsq[j], hlp);
                mexPrintf(temp);
#endif
                res[i+j*m] = exp(-(xsq[i]+xsq[j]-2*hlp)/gamma);
                /* Copy lower left triangular to upper right */
                res[j+i*m] = res[i+j*m];
            }
        }
        
    /* Two argument case */
    } else {
        

#ifdef DEBUG
    sprintf(temp, "Call with two parameters!\n");
    mexPrintf(temp);
#endif        
        if (n != mxGetM(prhs[2]))
            mexErrMsgTxt("GaussKernel.evaluate (mex): Array dimensions must match");

        k  = mxGetN(prhs[2]);
        y = mxGetPr(prhs[2]);
        
        plhs[0] = mxCreateDoubleMatrix(m, k, mxREAL);
        res = mxGetPr(plhs[0]);
        
        double ysq[k];
        sumsq(ysq, y, n, k);
        /* Check for the correct loop to split up, depending on whether more columns or
         rows are given */
        if (m > k) {
            #pragma omp parallel for shared(res,xsq,ysq,x,y,n,m,k) private(hlp,i,j,l)
            for (i = 0; i < m; i++) {
                for (j = 0; j < k; j++) {
                    hlp=0;
                    for (l = 0; l < n; l++) {
#ifdef DEBUG
                        sprintf(temp, "ScaP i=%d, j=%d, l=%d, x=%.12f, y=%.12f\n", i, j, l, x[i*n+l], y[j*n+l]);
                        mexPrintf(temp);
#endif
                        hlp += x[i*n+l]*y[j*n+l];
                    }
#ifdef DEBUG
                sprintf(temp, "xsq=%.12f, ysq=%.12f, x*y=%.12f\n", xsq[i], ysq[j], hlp);
                mexPrintf(temp);
#endif
                res[i+j*m] = exp(-(xsq[i]+ysq[j]-2*hlp)/gamma);
                }
            }
        } else {
            #pragma omp parallel for shared(res,xsq,ysq,x,y,n,m,k) private(hlp,i,j,l)
            for (j = 0; j < k; j++) {
                for (i = 0; i < m; i++) {
                    hlp=0;
                    for (l = 0; l < n; l++) {
#ifdef DEBUG
                        sprintf(temp, "ScaP i=%d, j=%d, l=%d, x=%.12f, y=%.12f\n", i, j, l, x[i*n+l], y[j*n+l]);
                        mexPrintf(temp);
#endif
                        hlp += x[i*n+l]*y[j*n+l];
                    }
#ifdef DEBUG
                sprintf(temp, "xsq=%.12f, ysq=%.12f, x*y=%.12f\n", xsq[i], ysq[j], hlp);
                mexPrintf(temp);
#endif
                res[i+j*m] = exp(-(xsq[i]+ysq[j]-2*hlp)/gamma);
                }
            }
        }
    }
}
