#include <mex.h>
#include <math.h>
#ifdef OMP
/*#include <omp.h>*/
#endif

#define Y(i,j) y[(j)*ys+i]

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    /*
     prhs[0] = ExplSolver object
     prhs[1] = ode function handle
     prhs[2] = times - row vector
     prhs[3] = x0 - column vector
     */
    double *times, *y, *x0, *tin, *yin, *yout, dt;
    
    int tidx, i, ys, tlen;
    
#ifdef DEBUG
    char temp[100];
#endif
    
    tlen  = mxGetN(prhs[2]);
    times = mxGetPr(prhs[2]);
    ys = mxGetM(prhs[3]);
    
#ifdef DEBUG
    printf("Initialized.\ntlen=%d, ysize=%d\n",tlen,ys);
#endif
    
    /* Create y result matrix */
    plhs[0] = mxCreateDoubleMatrix(ys, tlen, mxREAL);
    y = mxGetPr(plhs[0]);
#ifdef DEBUG
    printf("y/plhs[0] size: %dx%d\n",mxGetN(plhs[0]),mxGetM(plhs[0]));
#endif
    
    /* Write x0 to it*/
    x0 = mxGetPr(prhs[3]);
#ifdef OMP
    printf("Parallel switch on, num_threads=%d!\n",omp_get_num_threads());
    #pragma omp parallel for shared(y,x0) private(i)
#endif
    for(i=0;i<ys;i++) {
#ifdef DEBUG
        printf("Writing y(1,%d)=%.12f\n",i,x0[i]);
#endif
        y[i] = x0[i];
    }
    
    mxArray *odein[3];
    odein[0] = prhs[1];
    odein[1] = mxCreateDoubleMatrix(1, 1, mxREAL);
    odein[2] = mxCreateDoubleMatrix(ys, 1, mxREAL);
    tin = mxGetPr(odein[1]);
    yin = mxGetPr(odein[2]);
    
    mxArray *odeout[1];
    /*odeout[0] = mxCreateDoubleMatrix(ys, 1, mxREAL);*/
    
#ifdef DEBUG
    printf("Created CallMATLAB args\n");
#endif

    for(tidx = 1; tidx < tlen; tidx++) {
        /* Get time-step */
        dt = times[tidx]-times[tidx-1];
        /* Set ode inputs */
        tin[0] = times[tidx-1];
        
#ifdef DEBUG
        printf("tin=%.12f, dt=%.12f\n",times[tidx-1],dt);
#endif  
#ifdef OMP
        #pragma omp parallel for shared(yin,y,tidx) private(i)
#endif
        for (i=0;i<ys;i++) {
            yin[i] = Y(i,tidx-1);
#ifdef DEBUG
            printf("yin[%d]=%.12f\n",i,Y(i,tidx-1));
#endif
        }
        /* Evaluate ODE fun: odefun(times(idx-1),y(:,idx-1)) 
         int mexCallMATLAB(int nlhs, mxArray *plhs[], int nrhs,
  mxArray *prhs[], const char *functionName);
         */
#ifdef DEBUG
        printf("Calling matlab function\n");
#endif
        if (mexCallMATLAB(1, odeout, 3, odein, "feval") != 0)
            mexErrMsgTxt("Failed calling matlab function handle.");
        
        yout = mxGetPr(odeout[0]);
#ifdef DEBUG
        printf("Called matlab function\n");
#endif        
        /* Update vector */
#ifdef OMP
        #pragma omp parallel for shared(y,tidx,dt,yout) private(i)
#endif
        for (i=0;i<ys;i++) {
            Y(i,tidx) = Y(i,tidx-1) + dt*yout[i];
#ifdef DEBUG
            printf("Y(%d,%d)=%.12f, Y(i,tidx-1)=%.12f, yout[i]=%.12f\n",i,tidx,Y(i,tidx),Y(i,tidx-1),yout[i]);
#endif
        }
    }

#ifdef DEBUG
    printf("Destroying arrays\n");
#endif

    /*mxDestroyArray(odein[0]);
    mxDestroyArray(odein[1]);*/
}
