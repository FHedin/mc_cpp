#pragma OPENCL EXTENSION cl_khr_fp64 : enable

__kernel void NonBonded_full(__global const double *epsi, __global const double *sig, __global const double *q,
                             __global double *x, __global double *y, __global double *z,
                             __global int *expair, __global int *exlist,
                             __global double *elec, __global double *vdw, const uint natom)
{
  uint i = get_global_id(0); //global ID of this work item

  double e_buf = 0.;
  double v_buf = 0.;

  double ei = epsi[i];
  double qi = q[i];
  
  int k=0;
  int l=expair[i];
  
  for (uint j=(i+1); j<natom; j++)
  {
    
      if ((l>0) && k<l && (exlist[l+k] == j))
      {
        continue;
      }
    
      double dx = x[j] - x[i];
      double dy = y[j] - y[i];
      double dz = z[j] - z[i];

      double r2  = dx*dx + dy*dy + dz*dz;
      double r6  = r2*r2*r2;
      double r12 = r6*r6;
      double rt  = 1.0/sqrt(r2);

      double e = ei * epsi[j];
      double s = sig[i] + sig[j];
      
      e_buf += qi*q[j]*rt;
      v_buf += e*( (s/r12) - (s/r6) );

  }

  elec[i] = e_buf;
  vdw[i]  = v_buf;
}

