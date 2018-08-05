#include "mex.h"

void block_inverse(double* C, double* Ci, mwSize n)
{
   double a11, a12, a13, a22, a23, a33, det;
   int i;
   for(i = 0; i < n; i=i+9){        
        a11=C[i];a12=C[i+1];a13=C[i+2];a22=C[i+4];a23=C[i+5];a33=C[i+8];
        det=a11*a22*a33+a12*a23*a13+a13*a12*a23-a11*a23*a23-a12*a12*a33-a22*a13*a13;
        Ci[i]=(a22*a33-a23*a23)/det;
        Ci[i+1]=(a23*a13-a12*a33)/det;
        Ci[i+2]=(a12*a23-a13*a22)/det;
        Ci[i+3]=Ci[i+1];
        Ci[i+4]=(a11*a33-a13*a13)/det;
        Ci[i+5]=(a12*a13-a11*a23)/det;
        Ci[i+6]=Ci[i+2];
        Ci[i+7]=Ci[i+5];
        Ci[i+8]=(a11*a22-a12*a12)/det;        
   }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    if(nrhs!=1)
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","1 input required.");
    
    if(nlhs!=1)
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","1 output required.");
                
    double* C=mxGetPr(prhs[0]);
    mwSize n=mxGetM(prhs[0]);
    plhs[0]=mxCreateDoubleMatrix(n, 1, mxREAL); double* Ci=mxGetPr(plhs[0]);        
    block_inverse(C, Ci, n); 
};
