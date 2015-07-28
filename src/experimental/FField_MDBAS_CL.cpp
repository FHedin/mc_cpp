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

#include <fstream>
#include <sstream>

#include <chrono>

using namespace std;
using namespace cl;

FField_MDBAS_CL::FField_MDBAS_CL(AtomList& _at_List, PerConditions& _pbc, Ensemble& _ens,
                                 string _cutMode, double _ctoff, double _cton, double _dcut)
    : FField(_at_List, _pbc, _ens, _cutMode, _ctoff, _cton, _dcut)
{
    init_CL();
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
            //cout << "\t Max W-Item : " << dv.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>() << endl;
            cout << "\t Extensions : " << dv.getInfo<CL_DEVICE_EXTENSIONS>() << endl;
            cout << endl;
        }

    }

}

void FField_MDBAS_CL::init_CL()
{
    //get all platforms (drivers)
    vector<Platform> all_platforms;
    Platform::get(&all_platforms);

    if(all_platforms.size()==0) {
        cout<<" No OpenCL platform (i.e. drivers) found. Check OpenCL installation!\n";
        exit(-1);
    }

    // build a list of available gpus for computation
    for(Platform pl : all_platforms)
    {
        vector<Device> all_devices;
        pl.getDevices(CL_DEVICE_TYPE_GPU, &all_devices);

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
      cout << "No GPU found ! Check OpenCL installation!\n" << endl;
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
        //cout << "\t Max W-Item : " << dv.getInfo<CL_DEVICE_MAX_WORK_ITEM_SIZES>() << endl;
        cout << "\t Extensions : " << dv.getInfo<CL_DEVICE_EXTENSIONS>() << endl;
    }

    cout << "Restricting to one GPU for the moment, the first of the list ... " << endl;
    def_platform = gpu_platforms[0];
    def_gpu = gpu_devices[0];
    
    // add those gpus to the calculation context
    cl_context = Context(def_gpu);

    //TODO: put to arrays the content of cl files from ./kernels
//     ifstream t("kernel/NonBonded_full.cl");
//     stringstream buffer;
//     buffer << t.rdbuf();
//     NonBonded_full = buffer.str();
//     
//     ifstream t2("kernel/NonBonded_switch.cl");
//     stringstream buffer2;
//     buffer2 << t2.rdbuf();
//     NonBonded_switch = buffer2.str();
//     
//     cout << "Content of NonBonded_full : " << endl;
//     cout << NonBonded_full << endl;
//     
//     cout << "Content of NonBonded_switch : " << endl;
//     cout << NonBonded_switch << endl;
    
    // prepare the openCL source code which is written for the moment in FField_MDBAS_CL.hpp
    cl_sources.push_back( {NonBonded_full.c_str(),NonBonded_full.length()} );
    //cl_sources.push_back( {NonBonded_switch.c_str(),NonBonded_switch.length()} );

    // compile the sources for our context elements, if errors print them and exit
    cl_program = Program(cl_context,cl_sources);
    if(cl_program.build({def_gpu})!=CL_SUCCESS)
    {
        cout<<" Error while compiling OpenCL routines : " << endl;
//         for(Device dv : gpu_devices)
//         {
        cout<<" BUILD LOG for device " << def_gpu.getInfo<CL_DEVICE_NAME>() << " : " << cl_program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(def_gpu) << "\n";
//         }
        exit(1);
    }

    //create queue to which we will push commands
    cl_queue = CommandQueue(cl_context,def_gpu);

    // create memory buffers on the device for storing data
    const int nAtom = ens.getN();
    // kernel will read from the following but never write to
    x   = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    y   = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    z   = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    epsi= Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    sig = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    q   = Buffer(cl_context,CL_MEM_READ_ONLY,sizeof(double)*nAtom);
    
    //kernel will read and write (TODO: write only ?)
    ener = Buffer(cl_context,CL_MEM_READ_WRITE,sizeof(double)*nAtom);
    
    // we can transmit already the epsi sig and q values as they are constant 
    const vector<double>& m_q = at_List.getChargevect();
    const vector<double>& m_sigma = at_List.getSigmavect();
    const vector<double>& m_epsi = at_List.getEpsilonvect();
    cl_queue.enqueueWriteBuffer(epsi,CL_TRUE,0,sizeof(double)*nAtom,m_epsi.data());
    cl_queue.enqueueWriteBuffer(sig,CL_TRUE,0,sizeof(double)*nAtom,m_sigma.data());
    cl_queue.enqueueWriteBuffer(q,CL_TRUE,0,sizeof(double)*nAtom,m_q.data());

    kernel_full = Kernel(cl_program,"NonBonded_full");
    kernel_full.setArg(0,epsi);
    kernel_full.setArg(1,sig);
    kernel_full.setArg(2,q);
    kernel_full.setArg(3,x);
    kernel_full.setArg(4,y);
    kernel_full.setArg(5,z);
    kernel_full.setArg(6,ener);
    
//     //run the kernel
//     cl::Kernel kernel_add=cl::Kernel(program,"simple_add");
//     kernel_add.setArg(0,buffer_A);
//     kernel_add.setArg(1,buffer_B);
//     kernel_add.setArg(2,buffer_C);
//     queue.enqueueNDRangeKernel(kernel_add,cl::NullRange,cl::NDRange(10),cl::NullRange);
//     queue.finish();
//
//     int C[10];
//     //read result C from the device to array C
//     queue.enqueueReadBuffer(buffer_C,CL_TRUE,0,sizeof(int)*10,C);
//
//     cout<<" result: \n";
//     for(int i=0; i<10; i++) {
//         cout<<C[i]<<" ";
//     }
//     cout << endl;

}


void FField_MDBAS_CL::clean_CL()
{
}

double FField_MDBAS_CL::getE()
{
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

double FField_MDBAS_CL::getEtot()
{
    // electrostatic and vdw are performed together for minimising computations
    computeNonBonded_full();
    computeNonBonded14();

    //     cout << "Electrostatic energy (kcal/mol) : " << this->elec / CONSTANTS::kcaltoiu << endl;
    //     cout << "Van der Waals energy (kcal/mol) : " << this->vdw / CONSTANTS::kcaltoiu << endl;

    // all the components of internal energy
    if ( nBond > 0 )
        computeEbond();
    //     cout << "Bonds energy (kcal/mol) : " << this->bond / CONSTANTS::kcaltoiu << endl;

    if ( nAngle > 0 )
        computeEang();
    //     cout << "Angles energy (kcal/mol) : " << this->ang / CONSTANTS::kcaltoiu << endl;

    if ( nUb > 0 )
        computeEub();
    //     cout << "Urey Bradley energy (kcal/mol) : " << this->ub / CONSTANTS::kcaltoiu << endl;

    if ( nDihedral > 0 )
        computeEdihe();
    //     cout << "Dihedrals Energy (kcal/mol) : " << this->dihe / CONSTANTS::kcaltoiu << endl;

    if ( nImproper > 0 )
        computeEimpr();
    //     cout << "Impropers energy (kcal/mol) : " << this->impr / CONSTANTS::kcaltoiu << endl;

    /* --- Other types of energies here --- */
    /**/

    pot = elec + vdw + bond + ang + ub + dihe + impr;
    tot = pot + kin;

    //     cout << "Potential energy (kcal/mol) : " << this->pot / CONSTANTS::kcaltoiu << endl;
    //     cout << "Kinetic energy (kcal/mol) : " << this->kin / CONSTANTS::kcaltoiu << endl;
    //     cout << "Total energy (kcal/mol) : " << this->tot / CONSTANTS::kcaltoiu << endl;

    return tot;
}

double FField_MDBAS_CL::getEswitch()
{
    // electrostatic and vdw are performed together for minimising computations
    computeNonBonded_switch();
    computeNonBonded14_switch();

    //     cout << "Electrostatic energy (kcal/mol) : " << this->elec / CONSTANTS::kcaltoiu << endl;
    //     cout << "Van der Waals energy (kcal/mol) : " << this->vdw / CONSTANTS::kcaltoiu << endl;

    // all the components of internal energy
    if ( nBond > 0 )
        computeEbond();
    //     cout << "Bonds energy (kcal/mol) : " << this->bond / CONSTANTS::kcaltoiu << endl;

    if ( nAngle > 0 )
        computeEang();
    //     cout << "Angles energy (kcal/mol) : " << this->ang / CONSTANTS::kcaltoiu << endl;

    if ( nUb > 0 )
        computeEub();
    //     cout << "Urey Bradley energy (kcal/mol) : " << this->ub / CONSTANTS::kcaltoiu << endl;

    if ( nDihedral > 0 )
        computeEdihe();
    //     cout << "Dihedrals Energy (kcal/mol) : " << this->dihe / CONSTANTS::kcaltoiu << endl;

    if ( nImproper > 0 )
        computeEimpr();
    //     cout << "Impropers energy (kcal/mol) : " << this->impr / CONSTANTS::kcaltoiu << endl;

    /* --- Other types of energies here --- */
    /**/

    pot = elec + vdw + bond + ang + ub + dihe + impr;
    tot = pot + kin;

    //     cout << "Potential energy (kcal/mol) : " << this->pot / CONSTANTS::kcaltoiu << endl;
    //     cout << "Kinetic energy (kcal/mol) : " << this->kin / CONSTANTS::kcaltoiu << endl;
    //     cout << "Total energy (kcal/mol) : " << this->tot / CONSTANTS::kcaltoiu << endl;

    return tot;
}

void FField_MDBAS_CL::computeNonBonded_full()
{
    double lelec = 0.;
    double lvdw = 0.;
    double di[3], dj[3];
    double qi, qj;
    double epsi, epsj;
    double sigi, sigj;
    double rt;
    bool exclude;

    const int nAtom = ens.getN();

    const vector<int>& exclPair = excl->getExclPair();
    const vector<vector<int>>& exclList = excl->getExclList();

    // #ifdef _OPENMP
    //     #pragma omp parallel default(none) private(di,dj,qi,qj,epsi,epsj,sigi,sigj,exclude,rt) shared(exclPair,exclList) reduction(+:lelec,lvdw)
    //     {
    //         #pragma omp for schedule(dynamic) nowait
    // #endif
    for (int i = 0; i < nAtom - 1; i++)
    {
        at_List.getCoords(i,di);
        qi = at_List.getCharge(i);
        epsi = at_List.getEpsilon(i);
        sigi = at_List.getSigma(i);

        int k = 0;

        for (int j = i + 1; j < nAtom; j++)
        {
            exclude = false;
            if ((exclPair[i]>0) && (exclList[i][k] == j))
            {
                exclude = true;
                k++;

                if (k >= exclPair[i])
                    k = exclPair[i] - 1;
            }

            double pelec = 0.;
            double pvdw = 0.;
            if (!exclude)
            {
                at_List.getCoords(j,dj);
                qj = at_List.getCharge(j);
                epsj = at_List.getEpsilon(j);
                sigj = at_List.getSigma(j);

                rt = Tools::distance2(di, dj, pbc);
                rt = sqrt(rt);
                rt = 1. / rt;
                pelec = computeEelec(qi, qj, rt);
                pvdw  = computeEvdw(epsi, epsj, sigi, sigj, rt);

                lelec += pelec;
                lvdw  += pvdw;

            } // if not exclude
        } // inner loop
    } // outer loop

    // #ifdef _OPENMP
    //     }
    // #endif

    this->elec =  CONSTANTS::chgcharmm * CONSTANTS::kcaltoiu * lelec;
    this->vdw = 4.0 * lvdw;

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
    int i, j, k, l;
    double lelec = 0., pelec;
    double levdw = 0., pvdw;
    double r, r2, rt;
    double di[3], dj[3];
    double qi, qj;
    double epsi, epsj;
    double sigi, sigj;

    const int nAtom = ens.getN();
    const double ctoff2 = cutoff*cutoff;
    const double cton2 = cuton*cuton;
    const double switch2 = 1./(Tools::X3<double>(ctoff2-cton2));

    const vector<int>& neighPair = excl->getNeighPair();
    const vector<int>& neighOrder = excl->getNeighOrder();
    const vector<vector<int>>& neighList = excl->getNeighList();

    //     ofstream stdf;
    //     stdf.open("std.txt",ios_base::out);
    //     stdf.precision(12);


    // #ifdef _OPENMP
    //     #pragma omp parallel default(none) private(i,j,k,l,di,dj,qi,qj,r,r2,rt,epsi,epsj,sigi,sigj,pelec,pvdw) shared(neighPair,neighOrder,neighList) reduction(+:lelec,levdw)
    //     {
    //         #pragma omp for schedule(dynamic) nowait
    // #endif
    for ( l = 0; l < nAtom; l++ )
    {
        i=neighOrder[l];

        at_List.getCoords(i,di);
        qi = at_List.getCharge(i);
        epsi = at_List.getEpsilon(i);
        sigi = at_List.getSigma(i);

        for ( k = 0; k < neighPair[i]; k++ )
        {
            j = neighList[i][k];
            at_List.getCoords(j,dj);
            qj = at_List.getCharge(j);
            epsj = at_List.getEpsilon(j);
            sigj = at_List.getSigma(j);

            r2 = Tools::distance2(di, dj, pbc);

            if ( r2 <= ctoff2 )
            {
                r = sqrt(r2);
                rt = 1. / r;

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



            } // end if r2

            //             stdf << i << '\t' << j << '\t' << pelec << '\t' << pvdw << endl;

        }// end loop neighList
    }// end loop natom

    // #ifdef _OPENMP
    //     }
    // #endif

    this->elec = CONSTANTS::chgcharmm * CONSTANTS::kcaltoiu * lelec;
    this->vdw = 4.0 * levdw;

    //     stdf.close();
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
