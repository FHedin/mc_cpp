/*
 *  mc_cpp : A Molecular Monte Carlo simulations software.
 *  Copyright (C) 2013  Florent Hedin
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FFIELD_MDBAS_CL_H
#define FFIELD_MDBAS_CL_H

#ifdef OPENCL_EXPERIMENTAL

// #define __CL_ENABLE_EXCEPTIONS
#define CL_USE_DEPRECATED_OPENCL_1_1_APIS

#include <CL/cl.hpp>

#include "FField.hpp"

#define TO_KB 1024.0
#define TO_MB (TO_KB*TO_KB)

class FField_MDBAS_CL : public FField
{
public:
  FField_MDBAS_CL(AtomList& _at_List, PerConditions& _pbc, Ensemble& _ens,
                    std::string _cutMode="switch", double _ctoff=12.0, double _cton=10.0, double _dcut=2.0);
  ~FField_MDBAS_CL();

  static void list_CL_Devices_GPU();
  
  virtual double getE();

protected:
  
  uint nAtom;
  //static const size_t local_work_size = 128;
  
  std::string NonBonded_full;
  std::string NonBonded_switch;
  
  // lists of platforms and gpus used for calculations
  std::vector<cl::Platform> gpu_platforms;
  std::vector<cl::Device> gpu_devices;
  cl::Platform def_platform;
  cl::Device   def_gpu;

  // the context manages the gpus which are in use
  cl::Context cl_context;
  
  // this vector-like object will contain the openCL source code
  cl::Program::Sources cl_sources;
  // the program object will compile the sources and be used when calculating energies
  cl::Program cl_program;
  
  // queue in charge of dispatching the kernel work to compute units of a given gpu
  cl::CommandQueue cl_queue;
  
  // buffer objects will contain data on the gpus
  // coordinates + LJ params + charges
  cl::Buffer g_x,g_y,g_z,g_epsi,g_sig,g_q;
  // where the enrgy for each atom will be stored
  cl::Buffer g_elec,g_vdw;
  // local vector containing energy for each atom, filled from gpu
  std::vector<double> l_elec, l_vdw;
  
  cl::Kernel kernel_full;
  cl::Kernel kernel_switch;
  
  virtual void init_CL();
  virtual void clean_CL();
  
  virtual double getEtot();
  virtual double getEswitch();
  
  virtual void computeNonBonded_full();
  virtual void computeNonBonded14();
  
  virtual void computeNonBonded_switch();
  virtual void computeNonBonded14_switch();
  
  inline double computeEelec(const double qi, const double qj, const double rt);
  inline double computeEvdw(const double epsi, const double epsj, const double sigi,
                            const double sigj, const double rt);
  virtual void computeEbond();
  virtual void computeEang();
  virtual void computeEub();
  virtual void computeEdihe();
  virtual void computeEimpr();
  
  // stores in destination content of file specified by fname
  static void readKernel(const char fname[], std::string& destination);
};

#endif

#endif // FFIELD_MDBAS_VECT_H

