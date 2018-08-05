#include <stdlib.h>
#include <stdio.h> 
#include <string.h>

int pcg(const int n, double * b, double * xx)
{
  
  
  double * px = new double[n];
  double * qx = new double[n];
  double * zx = new double[n];
  double * wx = new double[n];
  
  double rho(0), rho_1(0), alpha(0), beta(0);
  
  memcpy(rx,b,sizeof(double)*n);
  dblas::set_zero(n,xx);
  
  double bnorm = dblas::norm2(n,b);
  double res = bnorm;
  
  pcg_iter_count = 0;
  
  while(true){
    
    if (pcg_iter_count >= max_pcg_iter) break;
 //   if ((res < pcg_rtol*bnorm) and (pcg_iter_count >= min_pcg_iter)) break;
    if ((res < pcg_rtol*bnorm) and (pcg_iter_count >= 0)) break;
    
    // apply preconditioner
    apply_preconditioner(rx,zx);
    //block_multiply<size_c>(n_cams,P,rx,zx); //z = solve(M, r);
    
    rho = dblas::dot(n,rx, zx);
    if (pcg_iter_count == 0)
      memcpy(px,zx,sizeof(double)*n);
    else
      {
	beta = rho / rho_1;
	for(int i = 0; i < n; ++i)
	  px[i] = zx[i] + beta*px[i];
      };
    
    eval_Axb(px,qx);
    alpha = rho / dblas::dot(n,px, qx);
    dblas::axpy(n,alpha,px,xx);
    
    if (pcg_iter_count % 21 == 0){
      eval_Axb(xx,wx);
      for(int i = 0; i < n ; ++i)
	 rx[i] = b[i] - wx[i];
    }else{
      dblas::axpy(n,-alpha,qx,rx);
    }
    
    rho_1 = rho;
    res = dblas::norm2(n,rx);
    ++pcg_iter_count;
  }
  std::cout<<"pcg converged at iteration "<<pcg_iter_count<<"\n";
  delete []px;
  delete []qx;
  delete []zx;
  delete []wx;
  return 0;
};
