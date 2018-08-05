#include "mex.h"
#include "iostream.h"
#include "fstream.h"
using namespace std;

void read_bundle_output(char* filename, int iteration, double* output); 
{
   ifstream file;
   file.open(filename);
   
   int int_in;
   double double_in;
   
   file>>int_in;file>>double_in;file>>double_in;file>>double_in;   
   int i;
   for(i=0;i<iteration;i++)
   {
        file>>int_in;
        output[i*8+5]=int_in;
        file>>double_in;
        output[i*8+6]=double_in;
        file>>int_in;
        output[i*8+7]=int_in;
        file>>double_in;
        output[i*8]=double_in;
        file>>double_in;
        output[i*8+1]=double_in;
        file>>double_in;
        output[i*8+2]=double_in;
        file>>double_in;
        output[i*8+3]=double_in;
        file>>double_in;
        output[i*8+4]=double_in;
   }
}


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    if(nrhs!=2)
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","2 inputs required.");
    
    if(nlhs!=1)
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","1 output required.");
                
    char* filename=mxGetChars(prhs[0]);
    mwSize iteration=mxGetScalar(prhs[1]);
    /* Col 1: f value   Col 2: g value   Col 3: h   Col 4: rho    Col 5: mu    Col 6: PCG flag   Col 7: PCG residual   Col 8: PCG iteration */
    plhs[0]=mxCreateDoubleMatrix(iteration*8, 1, mxREAL); 
    double* output=mxGetPr(plhs[0]);        
    read_bundle_output(filename, iteration, output); 
};
