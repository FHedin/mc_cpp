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
#include <chrono>

#include "cl.hpp"

using namespace std;

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

void FField_MDBAS_CL::init_CL()
{
  //get all platforms (drivers)
  vector<cl::Platform> all_platforms;
  cl::Platform::get(&all_platforms);
  if(all_platforms.size()==0){
    cout<<" No platforms found. Check OpenCL installation!\n";
    exit(1);
  }
  cl::Platform default_platform=all_platforms[0];
  cout << "Using platform: "<<default_platform.getInfo<CL_PLATFORM_NAME>()<<"\n";
  
  //get default device of the default platform
  vector<cl::Device> all_devices;
  default_platform.getDevices(CL_DEVICE_TYPE_ALL, &all_devices);
  if(all_devices.size()==0){
    cout<<" No devices found. Check OpenCL installation!\n";
    exit(1);
  }
  cl::Device default_device=all_devices[0];
  cout<< "Using device: "<<default_device.getInfo<CL_DEVICE_NAME>()<<"\n";
  
  
  cl::Context context({default_device});
  
  cl::Program::Sources sources;
  
  // kernel calculates for each element C=A+B
  string kernel_code=
  "   void kernel simple_add(global const int* A, global const int* B, global int* C){       "
  "       C[get_global_id(0)]=A[get_global_id(0)]+B[get_global_id(0)];                 "
  "   }                                                                               ";
  sources.push_back({kernel_code.c_str(),kernel_code.length()});
  
  cl::Program program(context,sources);
  if(program.build({default_device})!=CL_SUCCESS){
    cout<<" Error building: "<<program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(default_device)<<"\n";
    exit(1);
  }
  
  
  // create buffers on the device
  cl::Buffer buffer_A(context,CL_MEM_READ_WRITE,sizeof(int)*10);
  cl::Buffer buffer_B(context,CL_MEM_READ_WRITE,sizeof(int)*10);
  cl::Buffer buffer_C(context,CL_MEM_READ_WRITE,sizeof(int)*10);
  
  int A[] = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
  int B[] = {0, 1, 2, 0, 1, 2, 0, 1, 2, 0};
  
  //create queue to which we will push commands for the device.
  cl::CommandQueue queue(context,default_device);
  
  //write arrays A and B to the device
  queue.enqueueWriteBuffer(buffer_A,CL_TRUE,0,sizeof(int)*10,A);
  queue.enqueueWriteBuffer(buffer_B,CL_TRUE,0,sizeof(int)*10,B);
  
  
  //run the kernel
  cl::KernelFunctor simple_add(cl::Kernel(program,"simple_add"),queue,cl::NullRange,cl::NDRange(10),cl::NullRange);
  simple_add(buffer_A,buffer_B,buffer_C);
  
  //alternative way to run the kernel
  /*
   cl::Kernel kernel_add=cl::Kernel(program,"simple_add");
   kernel_add.setArg(0,buffer_A);
   kernel_add.setArg(1,buffer_B);
   kernel_add.setArg(2,buffer_C);
   queue.enqueueNDRangeKernel(kernel_add,cl::NullRange,cl::NDRange(10),cl::NullRange);
   queue.finish();*/
  
  int C[10];
  //read result C from the device to array C
  queue.enqueueReadBuffer(buffer_C,CL_TRUE,0,sizeof(int)*10,C);
  
  cout<<" result: \n";
  for(int i=0;i<10;i++){
    cout<<C[i]<<" ";
  }
  
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
    Vec4d potVDW(0.),potELEC(0.);
    Vec4d ep_i,ep_j,sig_i,sig_j,q_i,q_j;
    Vec4d xi,yi,zi,xj,yj,zj;
    Vec4d r12,r6,r2,rt;
    Vec4d dx,dy,dz;

    const Vec4d zeroes(0.0);
    const Vec4d inf(numeric_limits<double>::infinity());

    const size_t psize = 4;
    size_t remaining,end;

    const size_t nAtom = ens.getN();
    const vector<double>& x = at_List.getXvect();
    const vector<double>& y = at_List.getYvect();
    const vector<double>& z = at_List.getZvect();
    const vector<double>& sigma = at_List.getSigmavect();

    //copies instead of references because we will modify by putting at 0 some elements for disabling because of the exclusion list
    vector<double> q;
    vector<double> epsi;

    const vector<int>& exclPair = excl->getExclPair();
    const vector<vector<int>>& exclList = excl->getExclList();

    // #ifdef _OPENMP
    //     #pragma omp parallel default(none) private(ep_i,ep_j,sig_i,sig_j,q_i,q_j,xi,yi,zi,xj,yj,zj,r12,r6,r2,rt,dx,dy,dz,remaining,end,q,epsi) firstprivate(potVDW,potELEC) shared(x,y,z,sigma,exclPair,exclList)
    //     {
    //         #pragma omp for schedule(dynamic) nowait
    // #endif
    for(int i=0; i<(nAtom-1); i++)
    {
        //copy of charges and lj epsilons
        q = vector<double>(at_List.getChargevect());
        epsi = vector<double>(at_List.getEpsilonvect());

        xi = Vec4d(x[i]);
        yi = Vec4d(y[i]);
        zi = Vec4d(z[i]);

        sig_i = Vec4d(sigma[i]);
        ep_i  = Vec4d(epsi[i]);
        q_i = Vec4d(q[i]);

        remaining = (nAtom-(i+1))%psize;
        end = nAtom - remaining;

        int k=0;
        for (int j = i + 1; j < nAtom; j++)
        {
            if ((exclPair[i]>0) && (exclList[i][k] == j))
            {
                q[j]=0.;
                epsi[j]=0.;

                k++;

                if (k >= exclPair[i])
                    k = exclPair[i] - 1;
            }
        }

        for(int j=i+1; j<end; j+=psize)
        {
            //             q_j   = Vec4d(q[j],q[j+1],q[j+2],q[j+3]);
            //             ep_j  = Vec4d(epsi[j],epsi[j+1],epsi[j+2],epsi[j+3]);
            //             sig_j = Vec4d(sigma[j],sigma[j+1],sigma[j+2],sigma[j+3]);

            q_j.load(q.data()+j);
            ep_j.load(epsi.data()+j);
            sig_j.load(sigma.data()+j);

            sig_j += sig_i;
            ep_j  *= ep_i;

            //square the sigmas and scale epsilon by 4
            sig_j = square(sig_j);

            //             xj = Vec4d(x[j],x[j+1],x[j+2],x[j+3]);
            //             yj = Vec4d(y[j],y[j+1],y[j+2],y[j+3]);
            //             zj = Vec4d(z[j],z[j+1],z[j+2],z[j+3]);

            xj.load(x.data()+j);
            yj.load(y.data()+j);
            zj.load(z.data()+j);

            dx = xi - xj;
            dy = yi - yj;
            dz = zi - zj;
            pbc.applyPBC(dx,dy,dz);
            r2 = square(dx) + square(dy) + square(dz);

            //electrostatics
            //get 1/r for elec
            rt = sqrt(r2);
            rt = 1.0/rt;

            rt *= q_i * q_j;
            potELEC += rt;

            //van der waals
            //div sigma2 by r2 and keep result in r2
            r2 = sig_j / r2 ;

            r6 = pow_const(r2,3);
            r12 = square(r6);

            r12 -= r6;
            r12 *= ep_j;

            potVDW += r12;

        }// j loop

        //         remaining=0;
        if(remaining>0)
        {
            //r2 = Vec4d(numeric_limits<double>::infinity());
            dx = zeroes;
            dy = zeroes;
            dz = zeroes;
            sig_j = zeroes;
            ep_j = zeroes;
            q_j = zeroes;

            size_t j = end;
            for(size_t k=0; k<remaining; k++)
            {
                dx.insert(k,x[i]-x[j+k]);
                dy.insert(k,y[i]-y[j+k]);
                dz.insert(k,z[i]-z[j+k]);
                sig_j.insert( k , sigma[j+k] );
                ep_j.insert( k , epsi[j+k] );
                q_j.insert( k , q[j+k] );
            }

            pbc.applyPBC(dx,dy,dz);
            r2 = square(dx) + square(dy) + square(dz);
            r2 = select(r2==0.0,inf,r2);

            sig_j += sig_i;
            ep_j *= ep_i;

            //square the sigmas and scale epsilon by 4
            sig_j *= sig_j;

            //electrostatics
            //get 1/r for elec
            rt = sqrt(r2);
            rt = 1.0/rt;

            rt *= q_i * q_j;
            potELEC += rt;

            //van de waals
            //div sigma2 by r2 and keep result in r2
            r2 = sig_j / r2 ;

            r6 = pow_const(r2,3);
            r12 = square(r6);

            r12 -= r6;
            r12 *= ep_j;

            potVDW += r12;

        }// remaining j loop

    }// i loop

    // #ifdef _OPENMP
    //     #pragma omp critical
    //     {
    // #endif
    this->vdw  = 4.0 * horizontal_add(potVDW);
    this->elec = CONSTANTS::chgcharmm * CONSTANTS::kcaltoiu * horizontal_add(potELEC);
    // #ifdef _OPENMP
    //     }
    // #endif

    // #ifdef _OPENMP
    //     }// parallel section
    // #endif

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

    Vec4d potVDW(0.),potELEC(0.);
    Vec4d ep_i,ep_j,sig_i,sig_j,q_i,q_j;
    Vec4d xi,yi,zi,xj,yj,zj;
    Vec4d r12,r6,r2,rt;
    Vec4d dx,dy,dz;

    const Vec4d ones(1.0);
    const Vec4d zeroes(0.0);
    const Vec4d minf(-1.0*numeric_limits<double>::infinity());

    Vec4d switchFunc;
    Vec4db test1 , test2;

    const size_t psize = 4;
    size_t remaining,end;

    const size_t nAtom = ens.getN();
    const Vec4d ctoff2(cutoff*cutoff);
    const Vec4d cton2(cuton*cuton);
    const Vec4d switch2(1./pow_const(ctoff2-cton2,3));

    const vector<double>& x = at_List.getXvect();
    const vector<double>& y = at_List.getYvect();
    const vector<double>& z = at_List.getZvect();

    const vector<double>& q = at_List.getChargevect();
    const vector<double>& sigma = at_List.getSigmavect();
    const vector<double>& epsi = at_List.getEpsilonvect();

    const vector<int>& neighPair = excl->getNeighPair();
    const vector<int>& neighOrder = excl->getNeighOrder();
    const vector<vector<int>>& neighList = excl->getNeighList();

    //     ofstream vecf;
    //     vecf.open("vec.txt",ios_base::out);
    //     vecf.precision(12);

    // #ifdef _OPENMP
    //     #pragma omp parallel default(none) shared(x,y,z,q,epsi,sigma,neighPair,neighOrder,neighList) firstprivate(potVDW,potELEC) \
    //                          private(ep_i,ep_j,sig_i,sig_j,q_i,q_j,xi,yi,zi,xj,yj,zj,r12,r6,r2,rt,dx,dy,dz,switchFunc,test1,test2,remaining,end)
    //     {
    //         #pragma omp for schedule(dynamic) nowait
    // #endif
    for ( size_t k = 0; k < nAtom; k++ )
    {
        const size_t i = neighOrder[k];

        xi = Vec4d(x[i]);
        yi = Vec4d(y[i]);
        zi = Vec4d(z[i]);

        ep_i  = Vec4d(epsi[i]);
        sig_i = Vec4d(sigma[i]);
        q_i   = Vec4d(q[i]);

        remaining = neighPair[i] % psize;
        end = neighPair[i] - remaining;

        for ( size_t l = 0; l < end; l+=4 )
        {
            const size_t j0 = neighList[i][l];
            const size_t j1 = neighList[i][l+1];
            const size_t j2 = neighList[i][l+2];
            const size_t j3 = neighList[i][l+3];

            //             xj.load(x.data()+j0);
            //             yj.load(y.data()+j0);
            //             zj.load(z.data()+j0);

            xj = Vec4d(x[j0],x[j1],x[j2],x[j3]);
            yj = Vec4d(y[j0],y[j1],y[j2],y[j3]);
            zj = Vec4d(z[j0],z[j1],z[j2],z[j3]);

            dx = xi - xj;
            dy = yi - yj;
            dz = zi - zj;

            pbc.applyPBC(dx,dy,dz);
            r2 = square(dx) + square(dy) + square(dz);

            test1 = (r2 <= ctoff2);
            if (horizontal_or(test1))
            {

                //                 ep_j.load(epsi.data()+j0);
                //                 sig_j.load(sigma.data()+j0);
                //                 q_j.load(q.data()+j0);

                ep_j  = Vec4d(epsi[j0],epsi[j1],epsi[j2],epsi[j3]);
                sig_j = Vec4d(sigma[j0],sigma[j1],sigma[j2],sigma[j3]);
                q_j   = Vec4d(q[j0],q[j1],q[j2],q[j3]);

                //lorentz-berthelot rules
                sig_j += sig_i;
                ep_j  *= ep_i;

                //square the sigmas
                sig_j = square(sig_j);

                //electrostatics
                //get 1/r for elec
                rt  = sqrt(r2);
                rt  = 1.0/rt;
                rt *= q_i * q_j;

                //van der waals
                //div sigma2 by r2
                r6 = sig_j / r2 ;
                r6 = pow_const(r6,3);
                r12 = square(r6);
                r12 -= r6;
                r12 *= ep_j;

                switchFunc = ones;

                test2 = (r2 > cton2) && test1;
                if(horizontal_or(test2))
                {
                    const Vec4d switch1 = ctoff2-r2 ;
                    const Vec4d switch3 = pow_const(switch1,2)*(ctoff2 + 2.*r2 - 3.*cton2)*switch2;
                    switchFunc = select(test2,switch3,ones);
                }

                rt  *= switchFunc;
                r12 *= switchFunc;

                potELEC = if_add(test1,potELEC,rt);
                potVDW  = if_add(test1,potVDW,r12);



            } //ctoff2 if

            //             vecf << i << '\t' << j0 << '\t' << rt[0] << '\t' << r12[0] << endl;
            //             vecf << i << '\t' << j1 << '\t' << rt[1] << '\t' << r12[1] << endl;
            //             vecf << i << '\t' << j2 << '\t' << rt[2] << '\t' << r12[2] << endl;
            //             vecf << i << '\t' << j3 << '\t' << rt[3] << '\t' << r12[3] << endl;

        } // end l loop based on remaining/end

        //         remaining=0;
        if(remaining>0)
        {
            dx = zeroes;
            dy = zeroes;
            dz = zeroes;
            sig_j = zeroes;
            ep_j = zeroes;
            q_j = zeroes;
            r2 = zeroes;

            size_t j = neighList[i][end];
            for(size_t m=0; m<remaining; m++)
            {
                dx.insert(m, x[i]-x[j+m]);
                dy.insert(m, y[i]-y[j+m]);
                dz.insert(m, z[i]-z[j+m]);
                sig_j.insert(m, sigma[j+m] );
                ep_j.insert(m, epsi[j+m] );
                q_j.insert(m, q[j+m] );
            }

            pbc.applyPBC(dx,dy,dz);

            r2 = square(dx) + square(dy) + square(dz);
            // put to infinity where no interactions
            r2 = select( (ep_j==0.0) && (q_j==0.0),minf,r2);

            test1 = (abs(r2) <= ctoff2);
            if (horizontal_or(test1))
            {
                //lorentz-berthelot rules
                sig_j += sig_i;
                ep_j *= ep_i;

                //square the sigmas
                sig_j = square(sig_j);

                //electrostatics
                //get 1/r for elec
                rt = sqrt(r2);
                rt = 1.0/rt;
                rt *= q_i * q_j;

                //van der waals
                //div sigma2 by r2 and keep result in r2
                r6 = sig_j / r2 ;

                r6 = pow_const(r6,3);
                r12 = square(r6);

                r12 -= r6;
                r12 *= ep_j;

                switchFunc = ones;

                test2 = (r2 > cton2) && test1;
                if(horizontal_or(test2))
                {
                    const Vec4d switch1 = ctoff2-r2 ;
                    const Vec4d switch3 = pow_const(switch1,2)*(ctoff2 + 2.*r2 - 3.*cton2)*switch2 ;
                    switchFunc = select(test2,switch3,ones);
                }

                rt  *= switchFunc;
                r12 *= switchFunc;

                potELEC = if_add(test1,potELEC,rt);
                potVDW  = if_add(test1,potVDW,r12);



            } //ctoff2 if

            //             vecf << "Extra start" << endl;
            //             for(size_t m=0; m<remaining; m++)
            //               vecf << i << '\t' << j+m << '\t' << rt[m] << '\t' << r12[m] << endl;
            //             vecf << "Extra end" << endl;

        }// loop on remaining atoms

    } // end k loop on natoms

    //     vecf.close();

    // #ifdef _OPENMP
    //         #pragma omp critical
    //         {
    // #endif
    this->vdw  = 4.0 * horizontal_add(potVDW);
    this->elec = CONSTANTS::chgcharmm * CONSTANTS::kcaltoiu * horizontal_add(potELEC);
    // #ifdef _OPENMP
    //         }
    // #endif


    // #ifdef _OPENMP
    //     }
    // #endif

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
