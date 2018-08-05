#include "mex.h"
#include "math.h"

void calculate_one_residual(const double* c, const double* x, const double* o, double* f){
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

void calculate_residual(const double* state_x, mwSize n_obs, mwSize n_cams, const double* cidx, const double* pidx, const double* obs, const double* obs_mask_precon, mwSize obs_mask_precon_length, const double* obs_mask, mwSize obs_mask_length, double* f, double* f2, mwSize n_inputs)
{
    int iter;
    if(n_inputs==6){
        for(iter=0;iter<n_obs;iter++)
        {
            int camera=(int)cidx[iter];
            int point=(int)pidx[iter];
            calculate_one_residual(state_x+camera*9, state_x+n_cams*9+point*3, obs+iter*2, f+iter*2);
        }    
    }else if(n_inputs==8){
        int obs_mask_precon_counter=0;
        
        for(iter=0;iter<n_obs;iter++)
        {
            int camera=(int)cidx[iter];
            int point=(int)pidx[iter];
            calculate_one_residual(state_x+camera*9, state_x+n_cams*9+point*3, obs+iter*2, f+iter*2);
            if((obs_mask_precon_counter<obs_mask_precon_length)&&(iter==obs_mask_precon[obs_mask_precon_counter])){
                calculate_one_residual(state_x+camera*9, state_x+n_cams*9+point*3, obs+iter*2, f2+iter*2);
                obs_mask_precon_counter++;
            }            
        }   
    }else if(n_inputs==10){
        int obs_mask_precon_counter=0;
        int obs_mask_counter=0;
        for(iter=0;iter<n_obs;iter++)
        {
            int camera=(int)cidx[iter];
            int point=(int)pidx[iter];            
            if((obs_mask_counter<obs_mask_length)&&(iter==obs_mask[obs_mask_counter])){
                calculate_one_residual(state_x+camera*9, state_x+n_cams*9+point*3, obs+iter*2, f+iter*2);
                obs_mask_counter++;
            }  
            if((obs_mask_precon_counter<obs_mask_precon_length)&&(iter==obs_mask_precon[obs_mask_precon_counter])){
                calculate_one_residual(state_x+camera*9, state_x+n_cams*9+point*3, obs+iter*2, f2+iter*2);
                obs_mask_precon_counter++;
            }  
        }   
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{   
    if(nlhs!=2)
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","2 outputs required");           
    double* state_x=mxGetPr(prhs[0]);
    mwSize n_cams=mxGetScalar(prhs[1]);
    mwSize n_obs=mxGetScalar(prhs[2]);
    double* cidx=mxGetPr(prhs[3]);
    double* pidx=mxGetPr(prhs[4]);
    double* obs=mxGetPr(prhs[5]);
    double* obs_mask_precon;
    double* obs_mask;
    mwSize obs_mask_precon_length=0;
    mwSize obs_mask_length=0;
    
    if(nrhs==8) {
        obs_mask_precon=mxGetPr(prhs[6]);
        obs_mask_precon_length=mxGetScalar(prhs[7]);
    }else if(nrhs==10) {
        obs_mask_precon=mxGetPr(prhs[6]);
        obs_mask=mxGetPr(prhs[7]);
        obs_mask_length=mxGetScalar(prhs[8]);
        obs_mask_precon_length=mxGetScalar(prhs[9]);
    }else if(nrhs!=6) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","6, 8 or 10 inputs required.");
    }
    
    plhs[0]=mxCreateDoubleMatrix(2*n_obs, 1, mxREAL);
    double* f=mxGetPr(plhs[0]);    
    double* f2; 
    if((nrhs==8)||(nrhs==10)) {             
        plhs[1]=mxCreateDoubleMatrix(2*n_obs, 1, mxREAL); 
        f2=mxGetPr(plhs[1]);
    }else{
        plhs[1]=mxCreateDoubleMatrix(1, 1, mxREAL); 
        f2=mxGetPr(plhs[1]);
    }
    
    calculate_residual(state_x, n_obs, n_cams, cidx, pidx, obs, obs_mask_precon, obs_mask_precon_length, obs_mask, obs_mask_length, f, f2, nrhs);       
}

