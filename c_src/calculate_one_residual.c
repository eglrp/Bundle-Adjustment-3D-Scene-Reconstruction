#include "mex.h"

void calculate_one_residual(double* c, double* x, double* o, double* f){
    double o1 = o[0];
    double o2 = o[1];
    double c1 = c[0];
    double c2 = c[1];
    double c3 = c[2];
    double c4 = c[3];
    double c5 = c[4];
    double c6 = c[5];
    double c7 = c[6];
    double c8 = c[7];
    double c9 = c[8];
    double x1 = x[0];
    double x2 = x[1];
    double x3 = x[2];
    double   t4 = c2*c2;
    double   t5 = c3*c3;
    double   t6 = c1*c1;
    double   t7 = t4+t5+t6+1.0/4.503599627370496E15;
    double   t8 = 1/t7;
    double   t9 = sqrt(t7);
    double   t10 = cos(t9);
    double   t11 = t10-1.0;
    double   t12 = sin(t9);
    double   t13 = 1/sqrt(t7);
    double   t15 = t4*t8;
    double   t16 = c2*t12*t13;
    double   t17 = c1*c3*t11*t8;
    double   t19 = t5*t8;
    double   t20 = c3*t12*t13;
    double   t21 = c1*c2*t11*t8;
    double   t30 = c4*c7;
    double   t31 = t15+t19;
    double   t32 = t11*t31;
    double   t33 = t32+1.0;
    double   t34 = c7*t33*x1;
    double   t35 = t20+t21;
    double   t36 = c7*t35*x2;
    double   t37 = t16-t17;
    double   t38 = c7*t37*x3;
    double   t14 = t30+t34-t36+t38;
    double   t18 = t6*t8;
    double   t22 = c1*t12*t13;
    double   t29 = c2*c3*t11*t8;
    double   t44 = c5*c7;
    double   t45 = t18+t19;
    double   t46 = t11*t45;
    double   t47 = t46+1.0;
    double   t48 = c7*t47*x2;
    double   t49 = t20-t21;
    double   t50 = c7*t49*x1;
    double   t51 = t22+t29;
    double   t52 = c7*t51*x3;
    double   t23 = t44+t48+t50-t52;
    double   t24 = t15+t18;
    double   t25 = t11*t24;
    double   t26 = t25+1.0;
    double   t27 = t26*x3;
    double   t28 = t16+t17;
    double   t39 = t14*t14;
    double   t40 = t22-t29;
    double   t41 = t40*x2;
    double   t54 = t28*x1;
    double   t42 = c6+t27+t41-t54;
    double   t43 = 1/(t42*t42);
    double   t53 = t23*t23;
    double   t55 = t43*t53;
    double   t56 = t39*t43;
    double   t57 = t55+t56;
    double   t58 = c9*t57;
    double   t59 = c8+t58;
    double   t60 = t57*t59;
    double   t61 = t60+1.0;
    double   t62 = 1/t42;
    f[0] = -o1-t14*t61*t62;
    f[1] = -o2-t23*t61*t62;
}


/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double* c;              /* input scalar */
    double* x;               /* 1xN input matrix */
    double* o;                /* size of matrix */
    double* f;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=3) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Three inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }    
    c = mxGetPr(prhs[0]);
    x = mxGetPr(prhs[1]);
    o = mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(2, 1, mxREAL);
    f = mxGetPr(plhs[0]);

    calculate_one_residual(c, x, o, f);
}
