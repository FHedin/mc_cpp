#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void NonBonded_full(__global const double *epsi, __global const double *sig, __global const double *q,
                             __global double *x, __global double *y, __global double *z,
                             __global double *elec, __global double *vdw, const uint natom)
{
//   uint i = get_global_id(0);
//   uint j = get_global_id(1);
// 
//   if (j > i)
//   {
//       double dx = x[j] - x[i];
//       double dy = y[j] - y[i];
//       double dz = z[j] - z[i];
//       
//       double r2  = dx*dx + dy*dy + dz*dz;
//       double r6  = r2*r2*r2;
//       double r12 = r6*r6;
//       double rt  = 1.0/sqrt(r2);
//       
//       double e   = epsi[i]*epsi[j];
//       double s   = sig[i]+sig[j];
//       
//       elec[i] += q[i]*q[j]*rt;
//       vdw[i]  += e*( (s/r12) - (s/r6) ) ;
//   }

  uint i = get_global_id(0); //global ID of this work item

  double r,r2,r6,r12;

  double e_buf = 0.0;
  double v_buf = 0.0;

  for (uint j=(i+1); j<natom; j++)
  {
    double dx = x[j] - x[i];
    double dy = y[j] - y[i];
    double dz = z[j] - z[i];
    
    double r2  = dx*dx + dy*dy + dz*dz;
    double r6  = r2*r2*r2;
    double r12 = r6*r6;
    double rt  = 1.0/sqrt(r2);
    
    double e   = epsi[i]*epsi[j];
    double s   = sig[i]+sig[j];
    
    e_buf += q[i]*q[j]*rt;
    v_buf += e*( (s/r12) - (s/r6) ) ;
  }

  elec[i] = e_buf;
  vdw[i]  = v_buf;
  
//   uint i = get_global_id(0);
// 
//   elec[i] = i;
//   vdw[i]  = i*i;
}
