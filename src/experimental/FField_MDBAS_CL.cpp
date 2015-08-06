/*
 *  mc_cpp : A Molecular Monte Carlo simulations software.
 *  Copyright (C) 2013-2015  Florent Hedin
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

#include "FField_MDBAS_CL.hpp"

#ifdef OPENCL_EXPERIMENTAL

#include <cstdlib>
#include <cstdio>

#include <chrono>
#include <numeric>

#include "Constants.hpp"

using namespace std;
using namespace cl;

FField_MDBAS_CL::FField_MDBAS_CL(AtomList& _at_List, PerConditions& _pbc, Ensemble& _ens,
                                 string _cutMode, double _ctoff, double _cton, double _dcut)
    : FField(_at_List, _pbc, _ens, _cutMode, _ctoff, _cton, _dcut)
{
    nAtom = ens.getN();
    //after gpu calculation energy will be copied back to this vector
    l_elec = vector<double>(nAtom,0.0);
    l_vdw  = vector<double>(nAtom,0.0);
}

FField_MDBAS_CL::~FField_MDBAS_CL()
{
    clean_CL();
}

void FField_MDBAS_CL::list_CL_Devices_GPU()
{
    cout << "Listing all OpenCL capable GPUs available on this machine : " << endl;

    //get all platforms (drivers)
    vector<Platform> all_platforms;
    Platform::get(&all_platforms);

    if(all_platforms.size()==0) {
        cout<<" No platforms found. Check OpenCL installation!\n";
        return;
    }


    for(Platform pl : all_platforms)
    {
        cout << endl;
        cout << "Available platform : " << pl.getInfo<CL_PLATFORM_NAME>() << endl;
        cout << "\t Vendor     : " << pl.getInfo<CL_PLATFORM_VENDOR>() << endl;
        cout << "\t Profile    : " << pl.getInfo<CL_PLATFORM_PROFILE>() << endl;
        cout << "\t Version    : " << pl.getInfo<CL_PLATFORM_VERSION>() << endl;
        cout << "\t Extensions : " << pl.getInfo<CL_PLATFORM_EXTENSIONS>() << endl;

        vector<Device> all_devices;
        pl.getDevices(CL_DEVICE_TYPE_GPU, &all_devices);

        if(all_devices.size()==0)
        {
            cout << "No GPU device available for platform " << pl.getInfo<CL_PLATFORM_NAME>() << endl;
            continue;
        }

        cout << "Found " << all_devices.size() << " devices for platform " << pl.getInfo<CL_PLATFORM_NAME>() << " : " << endl;
        for(Device dv : all_devices)
        {
            cout << "\t Device : " << dv.getInfo<CL_DEVICE_NAME>() << endl;
            cout << "\t Vendor     : " << dv.getInfo<CL_DEVICE_VENDOR>() << endl;
            cout << "\t Profile    : " << dv.getInfo<CL_DEVICE_PROFILE>() << endl;
            cout << "\t Version    : " << dv.getInfo<CL_DEVICE_VERSION>() << endl;
            cout << "\t CL Version : " << dv.getInfo<CL_DEVICE_OPENCL_C_VERSION>() << endl;
            cout << "\t Driver     : " << dv.getInfo<CL_DRIVER_VERSION>() << endl;
            cout << "\t Glob. mem. : " << dv.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()/TO_MB << " MB" << endl;
            cout << "\t Max. alloc : " << dv.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>()/TO_MB << " MB (at once)" << endl;
            cout << "\t Cnst. mem. : " << dv.getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>()/TO_KB << " KB" << endl;
            cout << "\t Loc. mem.  : " << dv.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>()/TO_KB << " KB" << endl;
            cout << "\t Extensions : " << dv.getInfo<CL_DEVICE_EXTENSIONS>() << endl;
            cout << endl;
        }

    }

}

/**
 * TODO Avoid using C functions
 */
void FField_MDBAS_CL::readKernel(const char fname[], string& destination)
{
  
  FILE *f = NULL;
  f = fopen(fname,"rt");
  if (f == NULL)
  {
    cout << "Error while opening kernel file : " << fname << endl;
    exit(-1);
  }

  // Determine file size
  fseek(f, 0, SEEK_END);
  long size = ftell(f);
  
  char* where = (char*) malloc (sizeof(char)*size+1);
  
  rewind(f);
  fread(where, sizeof(char), size, f);

  where[size]='\0';
  
  destination = string(where,size+1);
  
  fclose(f);
  f = NULL;
  free(where);
}

void FField_MDBAS_CL::init_CL()
{
    //for storing error code
    cl_int ret;
    
    //get all platforms (drivers)
    vector<Platform> all_platforms;
    Platform::get(&all_platforms);

    if(all_platforms.size()==0) {
        cerr << " No OpenCL platform (i.e. drivers) found. Check OpenCL installation!\n";
        exit(-1);
    }

    // build a list of available gpus for computation
    for(Platform pl : all_platforms)
    {
        vector<Device> all_devices;
        ret = pl.getDevices(CL_DEVICE_TYPE_GPU, &all_devices);
//         if(ret != CL_SUCCESS)
//         {
//           cerr << "Error while calling getDevices(...) : error code = " << ret << endl;
//         }

        if(all_devices.size()==0)
        {
            continue;
        }
        else
        {
            gpu_platforms.push_back(pl);
            gpu_devices.insert(gpu_devices.end(),all_devices.begin(),all_devices.end());
        }

    }
    
    if(gpu_devices.size()==0)
    {
      cerr << "No GPU found ! Check OpenCL installation!\n" << endl;
      exit(-1);
    }

    // print info aout found gpus
    cout << "Found " << gpu_devices.size() << " available device(s) for GPU based calculations : " << endl;
    for(Device dv : gpu_devices)
    {
        cout << "\t Device     : " << dv.getInfo<CL_DEVICE_NAME>() << endl;
        cout << "\t Vendor     : " << dv.getInfo<CL_DEVICE_VENDOR>() << endl;
        cout << "\t Profile    : " << dv.getInfo<CL_DEVICE_PROFILE>() << endl;
        cout << "\t Version    : " << dv.getInfo<CL_DEVICE_VERSION>() << endl;
        cout << "\t CL Version : " << dv.getInfo<CL_DEVICE_OPENCL_C_VERSION>() << endl;
        cout << "\t Driver     : " << dv.getInfo<CL_DRIVER_VERSION>() << endl;
        cout << "\t Glob. mem. : " << dv.getInfo<CL_DEVICE_GLOBAL_MEM_SIZE>()/TO_MB << " MB" << endl;
        cout << "\t Max. alloc : " << dv.getInfo<CL_DEVICE_MAX_MEM_ALLOC_SIZE>()/TO_MB << " MB (at once)" << endl;
        cout << "\t Cnst. mem. : " << dv.getInfo<CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE>()/TO_KB << " KB" << endl;
        cout << "\t Loc. mem.  : " << dv.getInfo<CL_DEVICE_LOCAL_MEM_SIZE>()/TO_KB << " KB" << endl;
        cout << "\t Extensions : " << dv.getInfo<CL_DEVICE_EXTENSIONS>() << endl;
    }

    cout << "Restricting to one GPU for the moment, the first of the list ... " << endl;
    def_platform = gpu_platforms[0];
    def_gpu = gpu_devices[0];
    
    // add those gpus to the calculation context
    cl_context = Context(def_gpu);

    // read the CL kernels from external files
    readKernel("kernels/NonBonded_full.cl",NonBonded_full);
    readKernel("kernels/NonBonded_switch.cl",NonBonded_switch);
    
//     cout << "Content of NonBonded_full : " << endl;
//     cout << NonBonded_full << endl;
//     
//     cout << "Content of NonBonded_switch : " << endl;
//     cout << NonBonded_switch << endl;
    
    // prepare the openCL source code
    cl_sources.push_back( {NonBonded_full.c_str(),NonBonded_full.length()} );
    //cl_sources.push_back( {NonBonded_switch.c_str(),NonBonded_switch.length()} );

    // compile the sources for our context elements, if errors print them and exit
    cl_program = Program(cl_context,cl_sources);
    if(cl_program.build({def_gpu},"-cl-std=CL1.1")!=CL_SUCCESS)
    {
        cerr << " Error while compiling OpenCL routines : " << endl;
//         for(Device dv : gpu_devices)
//         {
        cerr << " BUILD LOG for device " << def_gpu.getInfo<CL_DEVICE_NAME>() << " : " << cl_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(def_gpu) << "\n";
//         }
        exit(1);
    }

    //create queue to which we will push commands
    cl_queue = CommandQueue(cl_context,def_gpu);

    const vector<double>& l_epsi = at_List.getEpsilonvect();
    const vector<double>& l_sigma = at_List.getSigmavect();
    const vector<double>& l_q = at_List.getChargevect();
    const vector<int>& l_exclPair = excl->getExclPair();
    
    // TODO : avoid this useless list reorganisation
    const vector<vector<int>>& exclList = excl->getExclList();
    vector<int> l_exclList;
    for(vector<int> sub : exclList)
      l_exclList.insert(l_exclList.end(),sub.begin(),sub.end());
    
    // create memory buffers on the device for storing data
    // kernel will read from the following but never write to
    g_x   = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    g_y   = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    g_z   = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    g_epsi= Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    g_sig = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    g_q   = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    g_expair = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(int)*nAtom);
    g_exlist = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(int)*nAtom*l_exclList.size());
    
    //kernel will write
    g_elec = Buffer(cl_context,CL_MEM_WRITE_ONLY,sizeof(double)*nAtom);
    g_vdw  = Buffer(cl_context,CL_MEM_WRITE_ONLY,sizeof(double)*nAtom);
    
    // we can transmit already the epsi sig and q values as they are constant 
    cl_queue.enqueueWriteBuffer(g_epsi,CL_TRUE,0,sizeof(double)*nAtom,l_epsi.data());
    cl_queue.enqueueWriteBuffer(g_sig,CL_TRUE,0,sizeof(double)*nAtom,l_sigma.data());
    cl_queue.enqueueWriteBuffer(g_q,CL_TRUE,0,sizeof(double)*nAtom,l_q.data());
    //and also the lists
    cl_queue.enqueueWriteBuffer(g_expair,CL_TRUE,0,sizeof(int)*nAtom,l_exclPair.data());
    cl_queue.enqueueWriteBuffer(g_exlist,CL_TRUE,0,sizeof(int)*nAtom,l_exclList.data());

    kernel_full = Kernel(cl_program,"NonBonded_full");
    
    ret = kernel_full.setArg(0,g_epsi);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 0 of kernel : error code : " << ret << endl;
    }
    
    ret = kernel_full.setArg(1,g_sig);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 1 of kernel : error code : " << ret << endl;
    }
    
    ret = kernel_full.setArg(2,g_q);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 2 of kernel : error code : " << ret << endl;
    }
    
    ret = kernel_full.setArg(3,g_x);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 3 of kernel : error code : " << ret << endl;
    }
    
    ret = kernel_full.setArg(4,g_y);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 4 of kernel : error code : " << ret << endl;
    }
    
    ret = kernel_full.setArg(5,g_z);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 5 of kernel : error code : " << ret << endl;
    }
    
    ret = kernel_full.setArg(6,g_expair);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 6 of kernel : error code : " << ret << endl;
    }
    
    ret = kernel_full.setArg(7,g_exlist);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 7 of kernel : error code : " << ret << endl;
    }
    
    ret = kernel_full.setArg(8,g_elec);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 8 of kernel : error code : " << ret << endl;
    }
    
    ret = kernel_full.setArg(9,g_vdw);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 9 of kernel : error code : " << ret << endl;
    }
    
    ret = kernel_full.setArg(10,sizeof(uint),&nAtom);
    if (ret != CL_SUCCESS)
    {
      cerr << "Error whith arg 10 of kernel : error code : " << ret << endl;
    }
    
    this->clInitialised=true;

}


void FField_MDBAS_CL::clean_CL()
{
}

double FField_MDBAS_CL::getE()
{
    if(!clInitialised)
      init_CL();
    
    double ener=0.0;

    cout << "Hi from getE of " << __FILE__ << endl;

    switch(this->cutMode)
    {
    case FULL:
        ener=getEtot();
        break;

    case SWITCH:
        ener=getEswitch();
        break;
    default:
        cerr << "Error : bad type of cutMode. file " << __FILE__ << " line " << __LINE__ << endl;
        exit(-100);
        break;
    }

    return ener;
}

void FField_MDBAS_CL::getE(double ener[10])
{

  if(!clInitialised)
      init_CL();
      
  cout << "Hi from getE of " << __FILE__ << endl;

    switch(this->cutMode)
    {
    case FULL:
        ener[0]=getEtot();
        break;

    case SWITCH:
        ener[0]=getEswitch();
        break;
    default:
        cerr << "Error : bad type of cutMode. file " << __FILE__ << " line " << __LINE__ << endl;
        exit(-100);
        break;
    }

    ener[1]=pot;
    ener[2]=kin;
    ener[3]=elec;
    ener[4]=vdw;
    ener[5]=bond;
    ener[6]=ang;
    ener[7]=ub;
    ener[8]=dihe;
    ener[9]=impr;
}

double FField_MDBAS_CL::getEtot()
{
    // electrostatic and vdw are performed together for minimising computations
    computeNonBonded_full();
    computeNonBonded14();

        cout << "Electrostatic energy (kcal/mol) : " << this->elec / CONSTANTS::kcaltoiu << endl;
        cout << "Van der Waals energy (kcal/mol) : " << this->vdw / CONSTANTS::kcaltoiu << endl;

    // all the components of internal energy
    if ( nBond > 0 )
        computeEbond();
        cout << "Bonds energy (kcal/mol) : " << this->bond / CONSTANTS::kcaltoiu << endl;

    if ( nAngle > 0 )
        computeEang();
        cout << "Angles energy (kcal/mol) : " << this->ang / CONSTANTS::kcaltoiu << endl;

    if ( nUb > 0 )
        computeEub();
        cout << "Urey Bradley energy (kcal/mol) : " << this->ub / CONSTANTS::kcaltoiu << endl;

    if ( nDihedral > 0 )
        computeEdihe();
        cout << "Dihedrals Energy (kcal/mol) : " << this->dihe / CONSTANTS::kcaltoiu << endl;

    if ( nImproper > 0 )
        computeEimpr();
        cout << "Impropers energy (kcal/mol) : " << this->impr / CONSTANTS::kcaltoiu << endl;

    /* --- Other types of energies here --- */
    /**/

    pot = elec + vdw + bond + ang + ub + dihe + impr;
    tot = pot + kin;

    cout << "Potential energy (kcal/mol) : " << this->pot / CONSTANTS::kcaltoiu << endl;
    cout << "Kinetic energy (kcal/mol) : " << this->kin / CONSTANTS::kcaltoiu << endl;
    cout << "Total energy (kcal/mol) : " << this->tot / CONSTANTS::kcaltoiu << endl;

    return tot;
}

double FField_MDBAS_CL::getEswitch()
{
    // electrostatic and vdw are performed together for minimising computations
    computeNonBonded_switch();
    computeNonBonded14_switch();

        cout << "Electrostatic energy (kcal/mol) : " << this->elec / CONSTANTS::kcaltoiu << endl;
        cout << "Van der Waals energy (kcal/mol) : " << this->vdw / CONSTANTS::kcaltoiu << endl;

    // all the components of internal energy
    if ( nBond > 0 )
        computeEbond();
        cout << "Bonds energy (kcal/mol) : " << this->bond / CONSTANTS::kcaltoiu << endl;

    if ( nAngle > 0 )
        computeEang();
        cout << "Angles energy (kcal/mol) : " << this->ang / CONSTANTS::kcaltoiu << endl;

    if ( nUb > 0 )
        computeEub();
        cout << "Urey Bradley energy (kcal/mol) : " << this->ub / CONSTANTS::kcaltoiu << endl;

    if ( nDihedral > 0 )
        computeEdihe();
        cout << "Dihedrals Energy (kcal/mol) : " << this->dihe / CONSTANTS::kcaltoiu << endl;

    if ( nImproper > 0 )
        computeEimpr();
        cout << "Impropers energy (kcal/mol) : " << this->impr / CONSTANTS::kcaltoiu << endl;

    /* --- Other types of energies here --- */
    /**/

    pot = elec + vdw + bond + ang + ub + dihe + impr;
    tot = pot + kin;

    cout << "Potential energy (kcal/mol) : " << this->pot / CONSTANTS::kcaltoiu << endl;
    cout << "Kinetic energy (kcal/mol) : " << this->kin / CONSTANTS::kcaltoiu << endl;
    cout << "Total energy (kcal/mol) : " << this->tot / CONSTANTS::kcaltoiu << endl;

    return tot;
}

void FField_MDBAS_CL::computeNonBonded_full()
{
    cl_int ret;
    cout << "Hi from computeNonBonded_full of " << __FILE__ << endl;
  
    const vector<double>& l_x = at_List.getXvect();
    const vector<double>& l_y = at_List.getYvect();
    const vector<double>& l_z = at_List.getZvect();
  
    cl_queue.enqueueWriteBuffer(g_x,CL_TRUE,0,sizeof(double)*nAtom,l_x.data());
    cl_queue.enqueueWriteBuffer(g_y,CL_TRUE,0,sizeof(double)*nAtom,l_y.data());
    cl_queue.enqueueWriteBuffer(g_z,CL_TRUE,0,sizeof(double)*nAtom,l_z.data());
    
//     ret=cl_queue.enqueueNDRangeKernel(kernel_full,cl::NullRange,cl::NDRange(nAtom));
    ret=cl_queue.enqueueNDRangeKernel(kernel_full,cl::NullRange,cl::NDRange(nAtom));
    if (ret != CL_SUCCESS)
    {
      cerr << "Error at enqueueNDRangeKernel : error code : " << ret << endl;
      exit(-1);
    }
    
    ret=cl_queue.finish();
    if (ret != CL_SUCCESS)
    {
      cerr << "Error at finish : error code : " << ret << endl;
      exit(-1);
    }
    
    cl_queue.enqueueReadBuffer(g_elec,CL_TRUE,0,sizeof(double)*nAtom,l_elec.data());
    cl_queue.enqueueReadBuffer(g_vdw,CL_TRUE ,0,sizeof(double)*nAtom,l_vdw.data());

//     for(int i=0; i< nAtom; i++)
//     {
//       cout << l_elec[i] << '\t' << l_vdw[i] << endl;
//     }
    
    double telec = std::accumulate(l_elec.begin(),l_elec.end(),0.0);
    double tvdw  = std::accumulate(l_vdw.begin(),l_vdw.end(),0.0);

    this->elec = CONSTANTS::chgcharmm * CONSTANTS::kcaltoiu * telec;
    this->vdw  = 4.0 * tvdw;

}

void FField_MDBAS_CL::computeNonBonded14()
{

    int i, j, k;
    double lelec = 0., pelec;
    double levdw = 0., pvdw;
    double r, r2, rt;
    double di[3], dj[3];
    double qi, qj;
    double epsi, epsj;
    double sigi, sigj;

    const int nPair14 = excl->getNPair14();

    const vector<int>& neighList14 = excl->getNeighList14();

    for ( k = 0; k < nPair14; k++ )
    {
        i = neighList14[2 * k];
        j = neighList14[2 * k + 1];

        at_List.getCoords(i,di);
        qi = at_List.getCharge(i);
        epsi = at_List.getEpsilon14(i);
        sigi = at_List.getSigma14(i);

        at_List.getCoords(j,dj);
        qj = at_List.getCharge(j);
        epsj = at_List.getEpsilon14(j);
        sigj = at_List.getSigma14(j);

        // 23 FLOP
        r2 = Tools::distance2(di, dj, pbc);

        // 5 FLOP
        r = sqrt(r2);
        rt = 1. / r;

        // 30 FLOP
        pelec = computeEelec(qi, qj, rt);
        pvdw = computeEvdw(epsi, epsj, sigi, sigj, rt);

        // 2 FLOP
        lelec += pelec;
        levdw += pvdw;
    }

    // 2 FLOP
    this->elec +=  CONSTANTS::chgcharmm * CONSTANTS::kcaltoiu * lelec;
    this->vdw  += 4.0 * levdw;

}

void FField_MDBAS_CL::computeNonBonded_switch()
{
//     cout << "Hi from computeNonBonded_switch of " << __FILE__ << endl;
//     
//     this->elec = CONSTANTS::chgcharmm * CONSTANTS::kcaltoiu * 0.0;
//     this->vdw = 4.0 * 0.0;

}

void FField_MDBAS_CL::computeNonBonded14_switch()
{
    int i, j, k;
    double lelec = 0., pelec;
    double levdw = 0., pvdw;
    double r, r2, rt;
    double di[3], dj[3];
    double qi, qj;
    double epsi, epsj;
    double sigi, sigj;

    const int nPair14 = excl->getNPair14();
    const vector<int>& neighList14 = excl->getNeighList14();

    const double ctoff2 = cutoff*cutoff;
    const double cton2 = cuton*cuton;
    const double switch2 = 1./(Tools::X3<double>(ctoff2-cton2));

    for ( k = 0; k < nPair14; k++ )
    {
        i = neighList14[2 * k];
        j = neighList14[2 * k + 1];

        at_List.getCoords(i,di);
        qi = at_List.getCharge(i);
        epsi = at_List.getEpsilon14(i);
        sigi = at_List.getSigma14(i);

        at_List.getCoords(j,dj);
        qj = at_List.getCharge(j);
        epsj = at_List.getEpsilon14(j);
        sigj = at_List.getSigma14(j);

        r2 = Tools::distance2(di, dj, pbc);

        r = sqrt(r2);
        rt = 1. / r;

        if ( r2 <= ctoff2 )
        {
            pelec = computeEelec(qi, qj, rt);
            pvdw = computeEvdw(epsi, epsj, sigi, sigj, rt);

            double switchFunc = 1.0;

            if ( r2 > cton2 )
            {
                double switch1 = ctoff2-r2;
                switchFunc = Tools::X2<double>(switch1)*(ctoff2 + 2.*r2 - 3.*cton2)*switch2;
            }

            pelec *= switchFunc;
            pvdw  *= switchFunc;

            lelec += pelec;
            levdw += pvdw;
        } // end if r2 ctoff2

    } // end for 1,4

    this->elec +=  CONSTANTS::chgcharmm * CONSTANTS::kcaltoiu * lelec;
    this->vdw  += 4.0 * levdw;
}

double FField_MDBAS_CL::computeEelec(const double qi, const double qj, const double rt)
{
    return qi * qj * rt;
}

double FField_MDBAS_CL::computeEvdw(const double epsi, const double epsj, const double sigi,
                                    const double sigj, const double rt)
{
    return epsi * epsj * (Tools::X12<double>((sigi + sigj) * rt) - Tools::X6<double>((sigi + sigj) * rt));
}

void FField_MDBAS_CL::computeEbond()
{
    int i, j, ll;
    int type;
    double di[3], dj[3];
    double r0, k;
    double d;
    double ebond = 0.0;

    for ( ll = 0; ll < nBond; ll++ )
    {
        i = bndList[ll].getAt1();
        j = bndList[ll].getAt2();

        at_List.getCoords(i,di);
        at_List.getCoords(j,dj);
        d = Tools::distance2(di, dj, pbc);
        d = sqrt(d);

        r0 = bndList[ll].getR0();
        k = bndList[ll].getK();
        type = bndList[ll].getType();

        switch ( type )
        {
        case BHARM:
            ebond += 0.5 * k * Tools::X2<double>(d - r0);
            break;

        case BMORSE:
        {
            double beta, morsea, morseb;
            beta = bndList[ll].getBeta();
            morsea = exp(-beta * (d - r0));
            morseb = Tools::X2<double>(morsea);
            ebond += k * (morseb - 2. * morsea) + k;
        }
        break;

        default:
            ebond += 0.5 * k * Tools::X2<double>(d - r0);
            break;
        }
    }
    this->bond = ebond;
}

void FField_MDBAS_CL::computeEang()
{
    int i, j, k, ll;
    double di[3], dj[3], dk[3], dab[3], dbc[3];
    double rab, rbc, /*rabt, rbct,*/ cost, /*sint,*/ theta;
    double kst, theta0;
    double eang = 0.0;

    //     const double dbl_epsilon = numeric_limits<double>::epsilon();

    for ( ll = 0; ll < nAngle; ll++ )
    {
        i = angList[ll].getAt1();
        j = angList[ll].getAt2();
        k = angList[ll].getAt3();
        kst = angList[ll].getK();
        theta0 = angList[ll].getTheta0();

        at_List.getCoords(i,di);
        at_List.getCoords(j,dj);
        at_List.getCoords(k,dk);

        rab = Tools::distance2(di, dj, pbc, dab);
        rab = sqrt(rab);
        //         rabt = 1. / rab;

        rbc = Tools::distance2(dk, dj, pbc, dbc);
        rbc = sqrt(rbc);
        //         rbct = 1. / rbc;

        cost = (dab[0] * dbc[0] + dab[1] * dbc[1] + dab[2] * dbc[2]) / (rab * rbc);
        //         sint = max(dbl_epsilon, sqrt(1.0 - (cost * cost)));
        theta = acos(cost);

        eang += 0.5 * kst * Tools::X2<double>(theta - theta0);
    }
    this->ang = eang;
}

void FField_MDBAS_CL::computeEub()
{
    int i, j, ll;
    double di[3], dj[3];
    double r0, k;
    double d;
    double ebond = 0.0;

    for ( ll = 0; ll < nUb; ll++ )
    {
        i = ubList[ll].getAt1();
        j = ubList[ll].getAt2();

        at_List.getCoords(i,di);
        at_List.getCoords(j,dj);
        d = Tools::distance2(di, dj, pbc);
        d = sqrt(d);

        r0 = ubList[ll].getR0();
        k = ubList[ll].getK();

        ebond += 0.5 * k * Tools::X2<double>(d - r0);
    }
    this->ub = ebond;
}

void FField_MDBAS_CL::computeEdihe()
{
    int i, j, k, l, ll;
    double di[3], dj[3], dk[3], dl[3];
    double dab[3], dbc[3], dcd[3];
    double pb[3], pc[3];
    double rbc, rpb, rpc, r2pb, r2pc;
    double pbpc, cosp, sinp, phi;
    double edihe = 0.;
    double kst, phi0, mult;
    int /*order,*/ type;

    const double twopi = CONSTANTS::PI;
    const double dbl_epsilon = numeric_limits<double>::epsilon();

    for ( ll = 0; ll < nDihedral; ll++ )
    {
        i = diheList[ll].getAt1();
        j = diheList[ll].getAt2();
        k = diheList[ll].getAt3();
        l = diheList[ll].getAt4();
        kst = diheList[ll].getK();
        phi0 = diheList[ll].getPhi0();
        mult = diheList[ll].getMult();
        //         order = diheList[ll].getOrder();
        type = diheList[ll].getType();

        at_List.getCoords(i,di);
        at_List.getCoords(j,dj);
        at_List.getCoords(k,dk);
        at_List.getCoords(l,dl);

        Tools::vec_substract(dj, di, dab);
        pbc.applyPBC(dab);

        rbc = sqrt(Tools::distance2(dk, dj, pbc, dbc));

        Tools::vec_substract(dl, dk, dcd);
        pbc.applyPBC(dcd);

        // construct first dihedral vector
        pb[0] = dab[1] * dbc[2] - dab[2] * dbc[1];
        pb[1] = dab[2] * dbc[0] - dab[0] * dbc[2];
        pb[2] = dab[0] * dbc[1] - dab[1] * dbc[0];
        r2pb = Tools::X2<double>(pb[0]) + Tools::X2<double>(pb[1]) + Tools::X2<double>(pb[2]);
        rpb = sqrt(r2pb);

        // construct second dihedral vector
        pc[0] = dbc[1] * dcd[2] - dbc[2] * dcd[1];
        pc[1] = dbc[2] * dcd[0] - dbc[0] * dcd[2];
        pc[2] = dbc[0] * dcd[1] - dbc[1] * dcd[0];
        r2pc = Tools::X2<double>(pc[0]) + Tools::X2<double>(pc[1]) + Tools::X2<double>(pc[2]);
        rpc = sqrt(r2pc);

        // determine dihedral angle
        pbpc = pb[0] * pc[0] + pb[1] * pc[1] + pb[2] * pc[2];
        cosp = pbpc / (rpb * rpc);
        sinp = (dbc[0]*(pc[1] * pb[2] - pc[2] * pb[1]) + dbc[1]*(pb[0] * pc[2] - pb[2] * pc[0]) +
                dbc[2]*(pc[0] * pb[1] - pc[1] * pb[0])) / (rpb * rpc * rbc);
        phi = atan2(sinp, cosp);

        // avoid singularity in sinp
        if ( sinp >= 0. )
        {
            sinp = max(dbl_epsilon, fabs(sinp));
        }
        else
        {
            sinp = -(max(dbl_epsilon, fabs(sinp)));
        }

        // calculate potential energy
        switch ( type )
        {
        case DCOS: // cosine dihedral
            edihe += kst * (1. + cos(mult * phi - phi0));
            break;

        case DHARM: // harmonic dihedral
            phi = phi - phi0;
            phi = phi - Tools::rint(phi / twopi) * twopi;
            edihe += 0.5 * kst * (phi * phi);
            break;

        default:
            edihe += kst * (1. + cos(mult * phi - phi0));
            break;
        }

    } // end of for loop on dihedrals
    this->dihe = edihe;
}

void FField_MDBAS_CL::computeEimpr()
{
    int i, j, k, l, ll;
    double di[3], dj[3], dk[3], dl[3];
    double dab[3], dbc[3], dcd[3];
    double pb[3], pc[3];
    double rbc, rpb, rpc, r2pb, r2pc;
    double pbpc, cosp, sinp, phi;
    double eimpr = 0.;
    double kst, phi0, mult;
    int /*order,*/ type;

    const double twopi = CONSTANTS::PI;
    const double dbl_epsilon = numeric_limits<double>::epsilon();

    for ( ll = 0; ll < nImproper; ll++ )
    {
        i = imprList[ll].getAt1();
        j = imprList[ll].getAt2();
        k = imprList[ll].getAt3();
        l = imprList[ll].getAt4();
        kst = imprList[ll].getK();
        phi0 = imprList[ll].getPhi0();
        mult = imprList[ll].getMult();
        //order = imprList[ll].getOrder();
        type = imprList[ll].getType();

        at_List.getCoords(i,di);
        at_List.getCoords(j,dj);
        at_List.getCoords(k,dk);
        at_List.getCoords(l,dl);

        Tools::vec_substract(dj, di, dab);
        pbc.applyPBC(dab);

        rbc = sqrt(Tools::distance2(dk, dj, pbc, dbc));

        Tools::vec_substract(dl, dk, dcd);
        pbc.applyPBC(dcd);

        // construct first dihedral vector
        pb[0] = dab[1] * dbc[2] - dab[2] * dbc[1];
        pb[1] = dab[2] * dbc[0] - dab[0] * dbc[2];
        pb[2] = dab[0] * dbc[1] - dab[1] * dbc[0];
        r2pb = Tools::X2<double>(pb[0]) + Tools::X2<double>(pb[1]) + Tools::X2<double>(pb[2]);
        rpb = sqrt(r2pb);

        // construct second dihedral vector
        pc[0] = dbc[1] * dcd[2] - dbc[2] * dcd[1];
        pc[1] = dbc[2] * dcd[0] - dbc[0] * dcd[2];
        pc[2] = dbc[0] * dcd[1] - dbc[1] * dcd[0];
        r2pc = Tools::X2<double>(pc[0]) + Tools::X2<double>(pc[1]) + Tools::X2<double>(pc[2]);
        rpc = sqrt(r2pc);

        // determine dihedral angle
        pbpc = pb[0] * pc[0] + pb[1] * pc[1] + pb[2] * pc[2];
        cosp = pbpc / (rpb * rpc);
        sinp = (dbc[0]*(pc[1] * pb[2] - pc[2] * pb[1]) + dbc[1]*(pb[0] * pc[2] - pb[2] * pc[0]) +
                dbc[2]*(pc[0] * pb[1] - pc[1] * pb[0])) / (rpb * rpc * rbc);
        phi = atan2(sinp, cosp);

        // avoid singularity in sinp
        if ( sinp >= 0. )
        {
            sinp = max(dbl_epsilon, fabs(sinp));
        }
        else
        {
            sinp = -(max(dbl_epsilon, fabs(sinp)));
        }

        // calculate potential energy
        switch ( type )
        {
        case DCOS: // cosine dihedral
            eimpr += kst * (1. + cos(mult * phi - phi0));
            break;

        case DHARM: // harmonic dihedral
            phi = phi - phi0;
            phi = phi - Tools::rint(phi / twopi) * twopi;
            eimpr += 0.5 * kst * (phi * phi);
            break;

        default:
            eimpr += kst * (1. + cos(mult * phi - phi0));
            break;
        }

    } // end of for loop on dihedrals
    this->impr = eimpr;
}

#endif /* VECTORCLASS_EXPERIMENTAL */
