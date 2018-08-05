#include "mex.h"
#include <stdio.h>
#include <string.h>

void set_zero(int n,  double *A)
  {
    /* to be optimized set the array to zero, use memset for better efficiency */
    int i;
    for(i = 0; i < n; ++i) A[i]=0;
  };
  
void mysymv(int n, double const *A, double const *b,double * c)
  {
    /* c = alpha A b+ beta c  */
    int i, j;
    for(i = 0; i<n; ++i){
    
      for(j = 0; j<i; ++j)
	c[i] += A[i*n+j]*b[j];
      
      for(j = i; j<n; ++j)
	c[i] += A[j*n+i]*b[j];
    };
  };

void  mygemvT(int n, int m, double const *A, double const *b,double * c)
  {
    /* c = c - A'b */
    int i, j;
    for(i = 0; i<m;++i){
      for(j = 0; j<n;++j){
	c[i]  -= A[i*n+j]*b[j];
      } 
    }
  }
    
void block_multiply(int n, int nblocks, double const *A, double const *b,double * c)
  {    
    set_zero(n*nblocks,c);
    int k;
    for (k = 0; k < nblocks; ++k){
      mysymv(n,A,b,c);
      A += n*n;
      b += n;
      c += n;
    }
  };
  
void build_b_large_scale(double* E, double* Ci, double* g, int n_cams, int n_obs, int n_pts, double* cidx, double* pidx, double* b)
{
    /* g=[p; q], tmpa=Ci*q */
    double* tmpa=mxGetPr(mxCreateDoubleMatrix(n_pts*3, 1, mxREAL));
    block_multiply(3, n_pts,Ci,g+9*n_cams,tmpa);    
    memcpy(b,g,sizeof(double)*n_cams*9);
    double* ptr = E; 
    /*  b = p - E' Ci q */
    int k;
    for(k = 0; k < n_obs; ++k){
      int cam = (int)cidx[k];
      int pt = (int)pidx[k];
      mygemvT(3, 9, ptr,tmpa +  pt*3, b + cam*9); 
      ptr += 27;
    };
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
        
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","8 inputs required.");
    }
    if(nlhs!=8) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","1 output required.");
    }    
    double* E=mxGetPr(prhs[0]);
    double* Ci=mxGetPr(prhs[1]);    
    double* g=mxGetPr(prhs[2]);
    mwSize n_cams=mxGetScalar(prhs[3]);
    mwSize n_obs=mxGetScalar(prhs[4]);
    mwSize n_pts=mxGetScalar(prhs[5]);
    double* cidx=mxGetPr(prhs[6]);
    double* pidx=mxGetPr(prhs[7]);    
    
    
    plhs[0]=mxCreateDoubleMatrix(81*n_cams, 1, mxREAL); 
    double* b=mxGetPr(plhs[0]);
        
    build_b_large_scale(E, Ci, g, n_cams, n_obs, n_pts, cidx, pidx, b);                   
};

