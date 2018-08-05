#include "mex.h"
#include "math.h"

void mysyrk(double const * A, double * B, int n, int rank)
{
/* in column major order, the upper triangle is non-zero */
    int i, j, k;
    for(i = 0; i< n;i++){
        for(j = i; j<n;j++){                               
	        for(k = 0;k<rank;k++)
	            B[j*n+i] += A[k*n+i]*A[k*n+j];
	    }       
	}
}

void mygemv(double const *A, double const *b,double * c, int n, int m)
{
    int i, j;
    for(i = 0;i<n;i++){
      for(j = 0; j<m;j++){
    c[i] -= A[j*n+i]*b[j];
      }
    }
}

void mygemmNT(double alpha, double const *A, double const *B,double * C, int n, int rank, int m)
  {
    int i, j, k;
    for (i = 0;i <n; i++){
        for(j = 0; j<m;j++){            
	        for (k =0;k<rank;k++)
	            C[j*n+i] += alpha*A[k*n+i]*B[k*m+j];
	    }        
	}           
  }
      
void calculate_one_jacobian(double* c, double* x, double* g)
{
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
    double   t64 = c2*c2;
    double   t65 = c3*c3;
    double   t66 = c1*c1;
    double   t67 = t64+t65+t66+1.0/4.503599627370496E15;
    double   t68 = 1/t67;
    double   t69 = sqrt(t67);
    double   t70 = cos(t69);
    double   t71 = t70-1.0;
    double   t72 = sin(t69);
    double   t73 = 1/sqrt(t67);
    double   t75 = t64*t68;
    double   t76 = c2*t72*t73;
    double   t77 = c1*c3*t68*t71;
    double   t79 = t65*t68;
    double   t80 = c3*t72*t73;
    double   t81 = c1*c2*t68*t71;
    double   t92 = c4*c7;
    double   t93 = t75+t79;
    double   t94 = t71*t93;
    double   t95 = t94+1.0;
    double   t96 = c7*t95*x1;
    double   t97 = t80+t81;
    double   t98 = c7*t97*x2;
    double   t99 = t76-t77;
    double   t100 = c7*t99*x3;
    double   t74 = t100+t92+t96-t98;
    double   t78 = t66*t68;
    double   t82 = c1*t72*t73;
    double   t89 = c2*c3*t68*t71;
    double   t123 = c5*c7;
    double   t124 = t78+t79;
    double   t125 = t124*t71;
    double   t126 = t125+1.0;
    double   t127 = c7*t126*x2;
    double   t128 = t80-t81;
    double   t129 = c7*t128*x1;
    double   t130 = t82+t89;
    double   t131 = c7*t130*x3;
    double   t83 = t123+t127+t129-t131;
    double   t84 = t75+t78;
    double   t85 = t71*t84;
    double   t86 = t85+1.0;
    double   t87 = t86*x3;
    double   t88 = t76+t77;
    double   t90 = 1/(t67*t67);
    double   t91 = 1/pow(t67,3.0/2.0);
    double   t101 = t74*t74;
    double   t102 = t82-t89;
    double   t103 = t102*x2;
    double   t133 = t88*x1;
    double   t104 = c6+t103-t133+t87;
    double   t105 = c1*t66*t90*2.0;
    double   t106 = c1*t64*t90*2.0;
    double   t139 = c1*t68*2.0;
    double   t107 = t105+t106-t139;
    double   t108 = t107*t71;
    double   t109 = c1*t72*t73*t84;
    double   t110 = t108+t109;
    double   t111 = c3*t66*t72*t91;
    double   t112 = c3*t66*t71*t90*2.0;
    double   t113 = c1*c2*t72*t91;
    double   t135 = c3*t68*t71;
    double   t136 = c1*c2*t68*t70;
    double   t114 = t111+t112+t113-t135-t136;
    double   t115 = t114*x1;
    double   t116 = t72*t73;
    double   t117 = t66*t68*t70;
    double   t118 = c1*c2*c3*t72*t91;
    double   t119 = c1*c2*c3*t71*t90*2.0;
    double   t138 = t66*t72*t91;
    double   t120 = t116+t117+t118+t119-t138;
    double   t121 = t120*x2;
    double   t145 = t110*x3;
    double   t122 = t115+t121-t145;
    double   t132 = t83*t83;
    double   t134 = 1/(t104*t104*t104);
    double   t137 = 1/(t104*t104);
    double   t140 = c1*t65*t90*2.0;
    double   t141 = c2*t66*t72*t91;
    double   t142 = c2*t66*t71*t90*2.0;
    double   t143 = c1*c3*t72*t91;
    double   t144 = t132*t137;
    double   t146 = t101*t122*t134*2.0;
    double   t147 = t122*t132*t134*2.0;
    double   t148 = t106+t140;
    double   t149 = t148*t71;
    double   t150 = c1*t72*t73*t93;
    double   t151 = t149+t150;
    double   t152 = c1*c3*t68*t70;
    double   t153 = t111+t112-t113-t135+t136;
    double   t154 = c7*t153*x3;
    double   t155 = -t116-t117+t118+t119+t138;
    double   t156 = c7*t155*x3;
    double   t157 = t105-t139+t140;
    double   t158 = t157*t71;
    double   t159 = c1*t124*t72*t73;
    double   t160 = t158+t159;
    double   t166 = c2*t68*t71;
    double   t161 = t141+t142-t143+t152-t166;
    double   t162 = c7*t161*x1;
    double   t274 = c7*t160*x2;
    double   t163 = t156+t162-t274;
    double   t164 = t101*t137;
    double   t165 = t144+t164;
    double   t167 = 1/t104;
    double   t168 = c9*t165;
    double   t169 = c8+t168;
    double   t170 = t165*t169;
    double   t171 = t170+1.0;
    double   t172 = c2*t64*t90*2.0;
    double   t173 = c2*t66*t90*2.0;
    double   t187 = c2*t68*2.0;
    double   t174 = t172+t173-t187;
    double   t175 = t174*t71;
    double   t176 = c2*t72*t73*t84;
    double   t177 = t175+t176;
    double   t178 = c3*t64*t72*t91;
    double   t179 = c3*t64*t71*t90*2.0;
    double   t180 = -t113-t135+t136+t178+t179;
    double   t181 = t180*x2;
    double   t182 = t64*t72*t91;
    double   t186 = t64*t68*t70;
    double   t183 = -t116+t118+t119+t182-t186;
    double   t184 = t183*x1;
    double   t193 = t177*x3;
    double   t185 = t181+t184-t193;
    double   t188 = c2*t65*t90*2.0;
    double   t189 = c1*t64*t72*t91;
    double   t190 = c1*t64*t71*t90*2.0;
    double   t191 = c2*c3*t68*t70;
    double   t192 = c2*c3*t72*t91;
    double   t194 = t173+t188;
    double   t195 = t194*t71;
    double   t196 = c2*t124*t72*t73;
    double   t197 = t195+t196;
    double   t210 = c1*t68*t71;
    double   t198 = t189+t190+t191-t192-t210;
    double   t199 = c7*t198*x1;
    double   t200 = t113-t135-t136+t178+t179;
    double   t201 = c7*t200*x3;
    double   t279 = c7*t197*x2;
    double   t202 = t199+t201-t279;
    double   t203 = t137*t202*t83*2.0;
    double   t204 = t116+t118+t119-t182+t186;
    double   t205 = c7*t204*x3;
    double   t206 = t172-t187+t188;
    double   t207 = t206*t71;
    double   t208 = c2*t72*t73*t93;
    double   t209 = t207+t208;
    double   t211 = t189+t190-t191+t192-t210;
    double   t212 = c7*t211*x2;
    double   t216 = c7*t209*x1;
    double   t213 = t205+t212-t216;
    double   t214 = t137*t213*t74*2.0;
    double   t277 = t101*t134*t185*2.0;
    double   t278 = t132*t134*t185*2.0;
    double   t215 = t203+t214-t277-t278;
    double   t217 = c1*t65*t72*t91;
    double   t218 = c1*t65*t71*t90*2.0;
    double   t219 = c3*t64*t90*2.0;
    double   t220 = -t191+t192-t210+t217+t218;
    double   t221 = t220*x1;
    double   t222 = c2*t65*t72*t91;
    double   t223 = c2*t65*t71*t90*2.0;
    double   t224 = -t143+t152-t166+t222+t223;
    double   t225 = t224*x2;
    double   t226 = c3*t66*t90*2.0;
    double   t227 = t219+t226;
    double   t228 = t227*t71;
    double   t229 = c3*t72*t73*t84;
    double   t230 = t228+t229;
    double   t242 = t230*x3;
    double   t231 = t221+t225-t242;
    double   t232 = t65*t72*t91;
    double   t233 = c3*t65*t90*2.0;
    double   t234 = t65*t68*t70;
    double   t241 = c3*t68*2.0;
    double   t235 = t219+t233-t241;
    double   t236 = t235*t71;
    double   t237 = c3*t72*t73*t93;
    double   t238 = t236+t237;
    double   t239 = t191-t192-t210+t217+t218;
    double   t240 = c7*t239*x3;
    double   t243 = t116+t118+t119-t232+t234;
    double   t244 = c7*t243*x1;
    double   t245 = t226+t233-t241;
    double   t246 = t245*t71;
    double   t247 = c3*t124*t72*t73;
    double   t248 = t246+t247;
    double   t249 = t143-t152-t166+t222+t223;
    double   t250 = c7*t249*x3;
    double   t283 = c7*t248*x2;
    double   t251 = t244+t250-t283;
    double   t252 = t137*t251*t83*2.0;
    double   t253 = -t116+t118+t119+t232-t234;
    double   t254 = c7*t253*x2;
    double   t286 = c7*t238*x1;
    double   t255 = t240+t254-t286;
    double   t256 = t137*t255*t74*2.0;
    double   t284 = t101*t134*t231*2.0;
    double   t285 = t132*t134*t231*2.0;
    double   t257 = t252+t256-t284-t285;
    double   t258 = t101*t134*2.0;
    double   t259 = t132*t134*2.0;
    double   t260 = t258+t259;
    double   t261 = t95*x1;
    double   t262 = t99*x3;
    double   t264 = t97*x2;
    double   t263 = c4+t261+t262-t264;
    double   t265 = t137*t263*t74*2.0;
    double   t266 = t126*x2;
    double   t267 = t128*x1;
    double   t299 = t130*x3;
    double   t268 = c5+t266+t267-t299;
    double   t269 = t137*t268*t83*2.0;
    double   t270 = t265+t269;
    double   t271 = t141+t142+t143-t152-t166;
    double   t272 = c7*t271*x2;
    double   t275 = c7*t151*x1;
    double   t273 = t154+t272-t275;
    double   t276 = t146+t147-t137*t163*t83*2.0-t137*t273*t74*2.0;
    double   t280 = t169*t215;
    double   t281 = c9*t165*t215;
    double   t282 = t280+t281;
    double   t287 = t169*t257;
    double   t288 = c9*t165*t257;
    double   t289 = t287+t288;
    double   t290 = c7*t137*t169*t74*2.0;
    double   t291 = c7*c9*t137*t165*t74*2.0;
    double   t292 = t290+t291;
    double   t293 = c7*t137*t169*t83*2.0;
    double   t294 = c7*c9*t137*t165*t83*2.0;
    double   t295 = t293+t294;
    double   t296 = t169*t260;
    double   t297 = c9*t165*t260;
    double   t298 = t296+t297;
    double   t300 = t169*t270;
    double   t301 = c9*t165*t270;
    double   t302 = t300+t301;
    double   t303 = t165*t165;
    double   t304 = t101*t134*t88*2.0;
    double   t305 = t132*t134*t88*2.0;
    double   t306 = c7*t128*t137*t83*2.0;
    double   t307 = c7*t137*t74*t95*2.0;
    double   t308 = t304+t305+t306+t307;
    double   t309 = t101*t102*t134*2.0;
    double   t310 = t102*t132*t134*2.0;
    double   t311 = c7*t137*t74*t97*2.0;
    double   t320 = c7*t126*t137*t83*2.0;
    double   t312 = t309+t310+t311-t320;
    double   t313 = t101*t134*t86*2.0;
    double   t314 = t132*t134*t86*2.0;
    double   t315 = c7*t130*t137*t83*2.0;
    double   t324 = c7*t137*t74*t99*2.0;
    double   t316 = t313+t314+t315-t324;
    double   t317 = t169*t308;
    double   t318 = c9*t165*t308;
    double   t319 = t317+t318;
    double   t321 = t169*t312;
    double   t322 = c9*t165*t312;
    double   t323 = t321+t322;
    double   t325 = t169*t316;
    double   t326 = c9*t165*t316;
    double   t327 = t325+t326;
    g[0] = -t167*t171*t273+t167*t74*((c8+c9*(t144+t101*1/pow(c6+t87-t88*x1+x2*(t82-c2*c3*t68*t71),2.0)))*(t146+t147-t137*t74*(t154-c7*t151*x1+c7*x2*(t141+t142+t143-c2*t68*t71-c1*c3*t68*t70))*2.0-t137*t163*t83*2.0)+c9*t165*(t146+t147-t137*t74*(t154+c7*x2*(t141+t142+t143-t152-c2*t68*t71)-c7*t151*x1)*2.0-t137*t163*t83*2.0))+t122*t137*t171*t74;
    g[1] = -t167*t171*t213-t167*t282*t74+t137*t171*t185*t74;
    g[2] = -t167*t171*(t240+c7*x2*(-t116+t118+t119+t232-t65*t68*t70)-c7*t238*x1)-t167*t289*t74+t137*t171*t231*t74;
    g[3] = -c7*t167*t171-t167*t292*t74;
    g[4] = -t167*t295*t74;
    g[5] = t137*t171*t74+t167*t298*t74;
    g[6] = -t167*t171*t263-t167*t302*t74;
    g[7] = -t165*t167*t74;
    g[8] = -t167*t303*t74;
    g[9] = t167*t83*(t169*t276+c9*t165*t276)-t163*t167*t171+t122*t137*t171*t83;
    g[10] = -t167*t171*t202-t167*t282*t83+t137*t171*t185*t83;
    g[11] = -t167*t171*t251-t167*t289*t83+t137*t171*t231*t83;
    g[12] = -t167*t292*t83;
    g[13] = -c7*t167*t171-t167*t295*t83;
    g[14] = t137*t171*t83+t167*t298*t83;
    g[15] = -t167*t171*t268-t167*t302*t83;
    g[16] = -t165*t167*t83;
    g[17] = -t167*t303*t83;
    g[18] = -t167*t319*t74-c7*t167*t171*t95-t137*t171*t74*t88;
    g[19] = t167*t323*t74+c7*t167*t171*t97+t102*t137*t171*t74;
    g[20] = t167*t327*t74-c7*t167*t171*t99+t137*t171*t74*t86;
    g[21] = -t167*t319*t83-c7*t128*t167*t171-t137*t171*t83*t88;
    g[22] = t167*t323*t83-c7*t126*t167*t171+t102*t137*t171*t83;
    g[23] = t167*t327*t83+c7*t130*t167*t171+t137*t171*t83*t86;
}

void calculate_jacobian(double * state_x, 
			       double * f,
			       int n_cams, int n_obs, int n_pts, double* cidx, double* pidx, double* D,
			       double * B, double * C, double * E, double * g)
{  
  double buffer[24];
  int iter;
  for(iter = 0; iter< n_obs; ++iter){
    int camera = (int)cidx[iter];
    int point = (int)pidx[iter];  

    calculate_one_jacobian(state_x+camera*9, state_x+n_cams*9+point*3, buffer);
   
    double * cweight = D + 9*camera;
    double * pweight= D + 9*n_cams + 3*point;
    int j, k;
    for (j = 0; j < 2; ++j){
	    for(k = 0; k < 9; ++k)
	        buffer[j*9+k] = buffer[j*9+k] * cweight[k];
	
	    for(k = 0; k < 3; ++k)
	        buffer[18 + j*3 + k]=buffer[18 + j*3 + k] *  pweight[k];
     }
        
    mysyrk(buffer,B+camera*81, 9, 2);
    mygemv(buffer,f+iter*2,g+camera*9, 9, 2);
    mygemmNT(1,buffer+18,buffer,E+iter*27, 3, 2, 9);
    mysyrk(buffer+18,C+point*9, 3, 2);
    mygemv(buffer+18,f+iter*2,g+n_cams*9+point*3, 3, 2);
  }
}

void calculate_camera_matrix(double * state_x, 
			       double * f,
			       int n_cams, int n_obs, int n_pts, double* cidx, double* pidx, double* D,
			       double* diag)
{
    double* B=new double[81*n_cams];
    double* C=new double[9*n_pts];
    double* Ci=new double[9*n_pts];
    double* E=new double[27*n_obs];
    double* g=new double[3*n_pts+9*n_cams];
    calculate_jacobian(state_x, f, n_cams, n_obs, n_pts, cidx, pidx, D, B, C, E, g);

    double min_diag=1e-6; double max_diag=1e+32;   
    double * ptr = B;
    for(int j=0;j<9;j++){
        double d = mu*std::min(std::max(ptr[j*10],min_diag),max_diag);   	            
        diag[j] = d;    
        ptr[j*10] += d; // add to A_ii
    }
    for(int i = 0; i < n_cams; ++i){
        for(int j = 0; j < 9; ++j){
            double d = mu*std::min(std::max(ptr[j*10],min_diag),max_diag);
            ptr[j*10] += d; // add to A_ii
        }
        ptr += 81;
    }       
    ptr = C; 
    for(int i = 0; i < n_pts; ++i){
        for(int j = 0; j < 3; ++j){
            double d = mu*std::min(std::max(ptr[j*4],min_diag),max_diag);      
            ptr[j*4] += d; // add to C_ii
        }
        ptr += 9;
    }

    // invert the point diagonal
    for(int i = 0; i < n_pts; ++i){
        a11=C[i];a12=C[i+3];a13=C[i+6];a22=C[i+4];a23=C[i+7];a33=C[i+8];
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
        ptr += 9;
    }
    
    for(int i=0;i<n_obs;i++){
        int camera = (int)cidx[i];
        int point = (int)pidx[i];  
        
    }
    delete B;
    delete C;
    delete Ci;
    delete E;
    delete g;    
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    if(nrhs!=8) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","8 inputs required.");
    }
    if(nlhs!=7) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","7 outputs required.");
    }    
    double* state_x=mxGetPr(prhs[0]);
    double* f=mxGetPr(prhs[1]);
    mwSize n_cams=mxGetScalar(prhs[2]);
    mwSize n_obs=mxGetScalar(prhs[3]);
    mwSize n_pts=mxGetScalar(prhs[4]);
    double* cidx=mxGetPr(prhs[5]);
    double* pidx=mxGetPr(prhs[6]);
    double* D=mxGetPr(prhs[7]);
    
    plhs[0]=mxCreateDoubleMatrix(81*n_cams, 1, mxREAL); 
    double* B=mxGetPr(plhs[0]);
    plhs[1]=mxCreateDoubleMatrix(9*n_pts, 1, mxREAL); 
    double* C=mxGetPr(plhs[1]);
    plhs[2]=mxCreateDoubleMatrix(27*n_obs, 1, mxREAL); 
    double* E=mxGetPr(plhs[2]);
    plhs[3]=mxCreateDoubleMatrix(3*n_pts+9*n_cams, 1, mxREAL); 
    double* g=mxGetPr(plhs[3]);
    plhs[4]=mxCreateDoubleMatrix(81*n_cams*2, 1, mxREAL);
    double* Bindices=mxGetPr(plhs[4]);
    plhs[5]=mxCreateDoubleMatrix(9*n_pts*2, 1, mxREAL);
    double* Cindices=mxGetPr(plhs[5]);
    plhs[6]=mxCreateDoubleMatrix(27*n_obs*2, 1, mxREAL);
    double* Eindices=mxGetPr(plhs[6]);
    
    calculate_jacobian(state_x, f, n_cams, n_obs, n_pts, cidx, pidx, D, B, C, E, g, Bindices, Cindices, Eindices);   
     
};

/*
#include <iostream>
#include <fstream>
#include <string>
using namespace std;
int main(int argc, char** argv ){
    int n_cams=49;
    int n_pts=7776;
    int n_obs=31843;
    double* state_x=new double[23769];
    double* f=new double[63686];
    int* cidx=new int[31843];
    int* pidx=new int[31843];
    double* D=new double[23769];
        
    ifstream input;
    input.open("state_x.txt");
    for(int i=0;i<23769;i++)
        input>>state_x[i];
    input.close();
    input.open("f.txt");
    for(int i=0;i<63686;i++)
        input>>f[i];
    input.close();
    input.open("cidx.txt");
    for(int i=0;i<31843;i++)
        input>>cidx[i];
    input.close();
   input.open("pidx.txt");
    for(int i=0;i<31843;i++)
        input>>pidx[i];
    input.close();
    input.open("D.txt");
    for(int i=0;i<23769;i++)
        input>>D[i];
    input.close();
    
    double* B=new double[3969];
    double* C=new double[69984];
    double* E=new double[859761];
    double* g=new double[23769];
    int* Bindices=new int[7938];
    int* Cindices=new int[139968];
    int* Eindices=new int[1719522];    
    calculate_jacobian(state_x, f, n_cams, n_obs, n_pts, cidx, pidx, D, B, C, E, g, Bindices, Cindices, Eindices); 
    return 0;
}
*/
