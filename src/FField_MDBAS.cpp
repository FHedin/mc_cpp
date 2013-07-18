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

#include <iomanip>
#include <algorithm> // for std::max
#include <limits> // for std::numeric_limits<double>
#include <chrono> // for precise timing
#include <string>

#include <cmath>
#include <cstdio>

// #include <papi.h>
// #include <scorep/SCOREP_User.h>

#include "FField_MDBAS.h"
#include "Constants.h"
#include "Tools.h"

using namespace std;

FField_MDBAS::FField_MDBAS(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens,
                           string _cutMode, double _ctoff, double _cton, double _dcut)
    : FField(_at_List, _pbc, _ens, _cutMode, _ctoff, _cton, _dcut)
{
}

FField_MDBAS::~FField_MDBAS()
{
}

double FField_MDBAS::getE()
{
    double ener=0.0;
    
//     float rtime;
//     float ptime;
//     long long flpops;
//     float mflops;
//     PAPI_flops(&rtime,&ptime,&flpops,&mflops);
    
//     long long ins;
//     float ipc;
//     PAPI_ipc(&rtime,&ptime,&ins,&ipc);
    
//     int Events[3] = {PAPI_L1_TCM,PAPI_L2_TCM,PAPI_L3_TCM};
//     long long int values[3] = {0,0,0};
//     PAPI_start_counters(Events,3);
    
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
    
//     PAPI_read_counters(values,3);
//     printf("L1 MISSES \t L2 MISSES \t L3 MISSES : %lld \t %lld \t %lld \n",values[0],values[1],values[2]);
//     PAPI_stop_counters(values,3);
    
//     PAPI_flops(&rtime,&ptime,&flpops,&mflops);
//     printf("MFLOP/s for FField_MDBAS::getE() : \t %f \n",mflops);

//     PAPI_ipc(&rtime,&ptime,&ins,&ipc);
//     printf("Average instructions per cycle for FField_MDBAS::getE() : \t %f \n",ipc);
    
    return ener;
}

double FField_MDBAS::getEtot()
{
//     cout << std::fixed << std::setprecision(15);

    // electrostatic and vdw are performed together for minimising computations
//     auto start = chrono::system_clock::now();
//     float rtime;
//     float ptime;
//     long long flpops;
//     float mflops;
//     PAPI_flops(&rtime,&ptime,&flpops,&mflops);
    
    computeNonBonded_full();
    computeNonBonded14();
    
//     PAPI_flops(&rtime,&ptime,&flpops,&mflops);
//     printf("Realtime\tCPUTime\tFLOP\tMFLOP/s\t%f\t%f\t%ld\t%f\n",rtime,ptime,flpops,mflops);
    
//     auto end = chrono::system_clock::now();
//     auto elapsed_time = chrono::duration_cast<chrono::nanoseconds> (end - start).count();
//     cout << "Electrostatic (kcal/mol) : " << this->elec / CONSTANTS::kcaltoiu << endl;
//     cout << "Van der Waals (kcal/mol) : " << this->vdw / CONSTANTS::kcaltoiu << endl;
//     cout << "Time required for NonBonded energy full was (nanoseconds) : " << elapsed_time << endl;
    
//     const int nAtom = ens.getN();
//     const int nPair14 = excl->getNPair14();
//     double flop_count = 60.0*(nAtom*(nAtom-1))/2.0;
//     flop_count += 60.0*nPair14 + 2.0;
//     cout << "Performances in GFLOPS for NonBonded energy full : " << (flop_count)/(double)(elapsed_time) << endl;


    // all the components of internal energy
//     start = chrono::system_clock::now();

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

//     end = chrono::system_clock::now();
//     elapsed_time = chrono::duration_cast<chrono::milliseconds> (end - start).count();
//     cout << "Time required for internal energy was (milliseconds) : " << elapsed_time << endl;

    /* --- Other types of energies here --- */

    pot = elec + vdw + bond + ang + ub + dihe + impr;
    tot = pot + kin;

//     cout << "Total energy (kcal/mol) : " << this->tot / CONSTANTS::kcaltoiu << endl;

    return tot;
}

double FField_MDBAS::getEswitch()
{
    // electrostatic and vdw are performed together for minimising computations
    computeNonBonded_switch();
    computeNonBonded14_switch();
    
    // all the components of internal energy
    if ( nBond > 0 )
        computeEbond();

    if ( nAngle > 0 )
        computeEang();

    if ( nUb > 0 )
        computeEub();

    if ( nDihedral > 0 )
        computeEdihe();

    if ( nImproper > 0 )
        computeEimpr();

    pot = elec + vdw + bond + ang + ub + dihe + impr;
    tot = pot + kin;

    return tot;
}

void FField_MDBAS::computeNonBonded_full()
{
//     SCOREP_USER_FUNC_BEGIN();
    
    int i, j, k;//, exclude=0;
    double lelec = 0.;//, pelec; 
    double levdw = 0.;//, pvdw;
//     double r, r2, rt;
//     double di[3], dj[3]; 
//     double qi, qj;
//     double epsi, epsj;
//     double sigi, sigj;

    // 4 BYTES copied
    const int nAtom = ens.getN();

    const vector<int>& exclPair = excl->getExclPair();
    const vector < vector<int >> &exclList = excl->getExclList();
    
    j=0;
    double* crds = new double[3*nAtom];
    double* q = new double[nAtom];
    double* e = new double[nAtom];
    double* s = new double[nAtom];
    
    for(i = 0; i < nAtom; i++)
    {
        at_List[i].getCoords(crds+j);
        q[i] = at_List[i].getCharge();
        e[i] = at_List[i].getEpsilon();
        s[i] = at_List[i].getSigma();
        j+=3;
    }

    for ( i = 0; i < nAtom - 1; i++ )
    {
        // 48 BYTES copied
//         at_List[i].getCoords(di); 
//         qi = at_List[i].getCharge();
//         epsi = at_List[i].getEpsilon();
//         sigi = at_List[i].getSigma();

        for ( j = i + 1; j < nAtom; j++ )
        {

            int exclude = 0;
            for ( k = 0; k < exclPair[i]; k++ )
            {
                if ( exclList[i][k] == j )
                {
                    exclude = 1;
                    break;
                }
            }

            if ( !exclude )
            {
                // 48 BYTES copied
//                 at_List[j].getCoords(dj);
//                 qj = at_List[j].getCharge();
//                 epsj = at_List[j].getEpsilon();
//                 sigj = at_List[j].getSigma();

                // 23 FLOP
                double r2 = Tools::distance2(crds+3*i, crds+3*j, pbc);

                // 4 FLOP (average with -O2 or -O3 optimisations)
                double r = sqrt(r2);
                
                // 1 FLOP
                double rt = 1. / r;

                // 4 FLOP
//                 pelec = computeEelec(qi, qj, rt);
                double pelec = computeEelec(q[i], q[j], rt);
                
                // 26 FLOP
//                 pvdw = computeEvdw(epsi, epsj, sigi, sigj, rt);
                double pvdw = computeEvdw(e[i], e[j], s[i], s[j], rt);
                
                // 2 FLOP
                lelec += pelec;
                levdw += pvdw;
            } // if not exclude
        } // inner loop
    } // outer loop

    this->elec = lelec;
    this->vdw = levdw;
    
    delete[] crds;
    delete[] q;
    delete[] e;
    delete[] s;
    
//     SCOREP_USER_FUNC_END();
}

void FField_MDBAS::computeNonBonded14()
{
//     SCOREP_USER_FUNC_BEGIN();
    
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

        at_List[i].getCoords(di);
        qi = at_List[i].getCharge();
        epsi = at_List[i].getEpsilon14();
        sigi = at_List[i].getSigma14();

        at_List[j].getCoords(dj);
        qj = at_List[j].getCharge();
        epsj = at_List[j].getEpsilon14();
        sigj = at_List[j].getSigma14();

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
    this->elec += lelec;
    this->vdw += levdw;
    
//     SCOREP_USER_FUNC_END();
}

// 4 FLOP
double FField_MDBAS::computeEelec(const double qi, const double qj, const double rt)
{
//     SCOREP_USER_FUNC_BEGIN();
    return CONSTANTS::chgcharmm * CONSTANTS::kcaltoiu * qi * qj * rt;
//     SCOREP_USER_FUNC_END();
}

// 26 FLOP
double FField_MDBAS::computeEvdw(const double epsi, const double epsj, const double sigi,
                                 const double sigj, const double rt)
{
//     SCOREP_USER_FUNC_BEGIN();
    return 4. * epsi * epsj * (Tools::X12<double>((sigi + sigj) * rt) - Tools::X6<double>((sigi + sigj) * rt));
//     SCOREP_USER_FUNC_END();
}

void FField_MDBAS::computeNonBonded_switch()
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

//     const vector<int>& exclPair = excl->getExclPair();
//     const vector<vector<int>>& exclList = excl->getExclList();

    const vector<int>& neighPair = excl->getNeighPair();
    const vector<int>& neighOrder = excl->getNeighOrder();
    const vector<vector<int>>& neighList = excl->getNeighList();

//     FILE *dataF=fopen("check.txt","wt");

    for ( l = 0; l < nAtom; l++ )
    {
        i=neighOrder[l];

        at_List[i].getCoords(di);
        qi = at_List[i].getCharge();
        epsi = at_List[i].getEpsilon();
        sigi = at_List[i].getSigma();

        for ( k = 0; k < neighPair[i]; k++ )
        {
            j = neighList[i][k];
            at_List[j].getCoords(dj);
            qj = at_List[j].getCharge();
            epsj = at_List[j].getEpsilon();
            sigj = at_List[j].getSigma();

            r2 = Tools::distance2(di, dj, pbc);

            r = sqrt(r2);
            rt = 1. / r;

//             fprintf(dataF,"%d\t%d\t%lf\t%lf\t%lf\n",i,j,r2,cton2,ctoff2);

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

            } // end if r2
        }// end loop neighList
    }// end loop natom

//     fclose(dataF);

    this->elec = lelec;
    this->vdw = levdw;
}

void FField_MDBAS::computeNonBonded14_switch()
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

        at_List[i].getCoords(di);
        qi = at_List[i].getCharge();
        epsi = at_List[i].getEpsilon14();
        sigi = at_List[i].getSigma14();

        at_List[j].getCoords(dj);
        qj = at_List[j].getCharge();
        epsj = at_List[j].getEpsilon14();
        sigj = at_List[j].getSigma14();

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

    this->elec += lelec;
    this->vdw += levdw;
}

#ifdef RANGED_E_EXPERIMENTAL

double FField_MDBAS::computeNonBonded_full_range(int first, int last)
{
    int i, j, k, exclude;
    double lelec = 0., pelec;
    double levdw = 0., pvdw;
    double r, r2, rt;
    double di[3], dj[3];
    double qi, qj;
    double epsi, epsj;
    double sigi, sigj;

    const int nAtom = ens.getN();

    const vector<int>& exclPair = excl->getExclPair();
    const vector < vector<int >> &exclList = excl->getExclList();

    for ( i = first; i <= last; i++ )
    {
        at_List[i].getCoords(di);
        qi = at_List[i].getCharge();
        epsi = at_List[i].getEpsilon();
        sigi = at_List[i].getSigma();

        for ( j = 0; j < nAtom; j++ )
        {
            if ( i == j )
                continue;

            exclude = 0;
            for ( k = 0; k < exclPair[i]; k++ )
            {
                if ( exclList[i][k] == j )
                {
                    exclude = 1;
                    break;
                }
            }

            if ( !exclude )
            {
                at_List[j].getCoords(dj);
                qj = at_List[j].getCharge();
                epsj = at_List[j].getEpsilon();
                sigj = at_List[j].getSigma();

                r2 = Tools::distance2(di, dj, pbc);

                r = sqrt(r2);
                rt = 1. / r;

                pelec = computeEelec(qi, qj, rt);
                pvdw = computeEvdw(epsi, epsj, sigi, sigj, rt);

                lelec += pelec;
                levdw += pvdw;
            }
        } // jloop
    } // i loop

    return (lelec + levdw);

}

double FField_MDBAS::computeNonBonded14_full_range(int first, int last)
{
    int i, j, k;
    double lelec = 0., pelec;
    double levdw = 0., pvdw;
    double r, r2, rt;
    double di[3], dj[3];
    double qi, qj;
    double epsi, epsj;
    double sigi, sigj;

    int nPair14 = excl->getNPair14();

    const vector<int>& neighList14 = excl->getNeighList14();

    for ( k = 0; k < nPair14; k++ )
    {
        i = neighList14[2 * k];
        j = neighList14[2 * k + 1];

//         cout << "in 1,4 ranged i,j,k  : " << i << '\t' << j << '\t' << k << '\t' << endl;

        if ( (i < first && j < first) || (i > last && j > last) )
            continue;

        at_List[i].getCoords(di);
        qi = at_List[i].getCharge();
        epsi = at_List[i].getEpsilon14();
        sigi = at_List[i].getSigma14();

        at_List[j].getCoords(dj);
        qj = at_List[j].getCharge();
        epsj = at_List[j].getEpsilon14();
        sigj = at_List[j].getSigma14();

        r2 = Tools::distance2(di, dj, pbc);

        r = sqrt(r2);
        rt = 1. / r;

        pelec = computeEelec(qi, qj, rt);
        pvdw = computeEvdw(epsi, epsj, sigi, sigj, rt);

        lelec += pelec;
        levdw += pvdw;
    }
    return (lelec + levdw);
}

#endif

void FField_MDBAS::computeEbond()
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

        at_List[i].getCoords(di);
        at_List[j].getCoords(dj);
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

void FField_MDBAS::computeEang()
{
    int i, j, k, ll;
    double di[3], dj[3], dk[3], dab[3], dbc[3];
    double rab, rbc, rabt, rbct, cost, sint, theta;
    double kst, theta0;
    double eang = 0.0;

    const double dbl_epsilon = numeric_limits<double>::epsilon();

    for ( ll = 0; ll < nAngle; ll++ )
    {
        i = angList[ll].getAt1();
        j = angList[ll].getAt2();
        k = angList[ll].getAt3();
        kst = angList[ll].getK();
        theta0 = angList[ll].getTheta0();

        at_List[i].getCoords(di);
        at_List[j].getCoords(dj);
        at_List[k].getCoords(dk);

        rab = Tools::distance2(di, dj, pbc, dab);
        rab = sqrt(rab);
        rabt = 1. / rab;

        rbc = Tools::distance2(dk, dj, pbc, dbc);
        rbc = sqrt(rbc);
        rbct = 1. / rbc;

        cost = (dab[0] * dbc[0] + dab[1] * dbc[1] + dab[2] * dbc[2]) / (rab * rbc);
        sint = max(dbl_epsilon, sqrt(1.0 - (cost * cost)));
        theta = acos(cost);

        eang += 0.5 * kst * Tools::X2<double>(theta - theta0);
    }
    this->ang = eang;
}

void FField_MDBAS::computeEub()
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

        at_List[i].getCoords(di);
        at_List[j].getCoords(dj);
        d = Tools::distance2(di, dj, pbc);
        d = sqrt(d);

        r0 = ubList[ll].getR0();
        k = ubList[ll].getK();

        ebond += 0.5 * k * Tools::X2<double>(d - r0);
    }
    this->ub = ebond;
}

void FField_MDBAS::computeEdihe()
{
    int i, j, k, l, ll;
    double di[3], dj[3], dk[3], dl[3];
    double dab[3], dbc[3], dcd[3];
    double pb[3], pc[3];
    double rbc, rpb, rpc, r2pb, r2pc;
    double pbpc, cosp, sinp, phi;
    double edihe = 0.;
    double kst, phi0, mult;
    int order, type;

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
        order = diheList[ll].getOrder();
        type = diheList[ll].getType();

        at_List[i].getCoords(di);
        at_List[j].getCoords(dj);
        at_List[k].getCoords(dk);
        at_List[l].getCoords(dl);

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

void FField_MDBAS::computeEimpr()
{
    int i, j, k, l, ll;
    double di[3], dj[3], dk[3], dl[3];
    double dab[3], dbc[3], dcd[3];
    double pb[3], pc[3];
    double rbc, rpb, rpc, r2pb, r2pc;
    double pbpc, cosp, sinp, phi;
    double eimpr = 0.;
    double kst, phi0, mult;
    int order, type;

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
        order = imprList[ll].getOrder();
        type = imprList[ll].getType();

        at_List[i].getCoords(di);
        at_List[j].getCoords(dj);
        at_List[k].getCoords(dk);
        at_List[l].getCoords(dl);

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

double FField_MDBAS::E_moving_set(int moveAtomList[], int moveBondList[])
{
    double ener = 0.0;

    // first get nonbonded energy for all moving atoms
    int ng = moveAtomList[0];
    int endng = ng + 2;
    int nn;
    int iaf, ial;

    for ( int it1 = 1; it1 <= ng; it1++ )
    {
        nn = moveAtomList[it1];
        for ( int it2 = endng; it2 <= nn; it2 += 2 )
        {
            iaf = moveAtomList[it2 - 1];
            ial = moveAtomList[it2];

//             cout << "iaf : " << iaf << " ial : " << ial << endl;
            
//             ener += ff.computeNonBonded_full_range(iaf, ial);
//             ener += ff.computeNonBonded14_full_range(iaf, ial);
        }
        endng = nn + 2;
    }

    return ener;
}