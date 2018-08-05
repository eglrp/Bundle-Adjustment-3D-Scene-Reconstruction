#include "mex.h"
#include "math.h"
#include <stdio.h>
#include <string.h>

void set_zero(int n, double* A)
{
  /*  for(int i = 0; i < n; ++i) A[i]=0; */
    memset(A, 0, n*sizeof(double));
}

void mygemvp(double const *A, double const *b,double * c, int n, int m)
  {
    /* c = c + Ab */
    int i, j;
    for(i = 0;i<n;++i){
      for(j = 0; j<m;++j){
    c[i] += A[j*n+i]*b[j];
      };
    };
  }


void mysymv(double const *A, double const *b,double * c, int n)
{
    /* c = alpha A b+ beta c  */
    int i, j;
    for(i = 0; i<n; ++i){
        for(j = 0; j<i; ++j)
            c[i] += A[i*n+j]*b[j];

        for(j = i; j<n; ++j)
            c[i] += A[j*n+i]*b[j];
    }
}

void mygemvT(double const *A, double const *b,double * c, int n, int m)
  {
    /* c = c - A'b */
    int i, j;
    for(i = 0; i<m;++i){
      for(j = 0; j<n;++j){
	c[i]  -= A[i*n+j]*b[j];
      } 
    }
  }
    
void block_multiply(int nblocks, double const *A, double const *b,double * c, int n)
{
/* multiply a vector with a block diagonal matrix */
    set_zero(n*nblocks,c);
    int k;
    for (k = 0; k < nblocks; ++k){
        mysymv(A,b,c,n);
        A += n*n;
        b += n;
        c += n;
    }
}

void eval_Axb(double *x, double *y, double* A, double* B, double* Ci, int n_cams, int n_pts, int n_obs, double* cidx, double* pidx, double* Bx, double* CiBx)
{
/*   y = (A - B'*Ci*B)x; */
    
    
    block_multiply(n_cams,A,x,y, 9);    
    
    set_zero(3*n_pts,Bx);
    double * ptr = B;
    int k;
    for(k = 0; k < n_obs; ++k){
        int cam = (int)cidx[k];
        int pt = (int)pidx[k];
        mygemvp(ptr, x +  cam*9, Bx + pt*3, 3, 9);
        ptr += 27;
    }
   
    block_multiply(n_pts,Ci,Bx,CiBx, 3);
    ptr = B; 
    for(k = 0; k < n_obs; ++k){
        int cam = cidx[k];
        int pt = pidx[k];
        mygemvT(ptr,CiBx +  pt*3, y + cam*9, 3, 9);
        ptr += 27;
    }  
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    if(nrhs!=9) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","9 inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","1 output required.");
    }
    double* x=mxGetPr(prhs[0]);
    double* A=mxGetPr(prhs[1]);
    double* B=mxGetPr(prhs[2]);
    double* Ci=mxGetPr(prhs[3]);
    mwSize n_cams=mxGetScalar(prhs[4]);
    mwSize n_pts=mxGetScalar(prhs[5]);
    mwSize n_obs=mxGetScalar(prhs[6]);
    double* cidx=mxGetPr(prhs[7]);
    double* pidx=mxGetPr(prhs[8]);     
    
    plhs[0]=mxCreateDoubleMatrix(9*n_cams, 1, mxREAL); 
    double* y=mxGetPr(plhs[0]);    
    
    double* Bx = mxGetPr(mxCreateDoubleMatrix(3*n_pts, 1, mxREAL));       
    double* CiBx = mxGetPr(mxCreateDoubleMatrix(3*n_pts, 1, mxREAL));
    
    eval_Axb(x, y, A, B, Ci, n_cams, n_pts, n_obs, cidx, pidx, Bx, CiBx);     
}
