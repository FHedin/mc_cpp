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

#include <cstdio>
#include <iostream>
#include <tuple>

#include "MC_metropolis.h"
#include "Move_TRN.h"
#include "Move_ROT.h"

using namespace std;

MC_metropolis::MC_metropolis(std::vector<Atom>& _at_List, PerConditions& _pbc,
                             Ensemble& _ens, FField& _ff, List_Moves& _mvlist, int _steps, double _dmax,
                             int _save_freq) : MC(_at_List, _pbc, _ens, _ff, _mvlist)
{
    nsteps = _steps;
    dmax = _dmax;
    svFreq = _save_freq;

    xyz = nullptr;
    xyz = fopen("tr.xyz", "w");

    efile = nullptr;
    efile = fopen("ener.dat", "w");

    cout << "Initialising MC Metropolis simulation : found " << ens.getN() << " atoms. The ensemble is " << ens.whoami() << std::endl;
}

MC_metropolis::~MC_metropolis()
{
    fclose(xyz);
    fclose(efile);
}

void MC_metropolis::run()
{
    int nmvtyp = mvlist.getNMoveTypes();
    if ( nmvtyp == 0 )
    {
        cerr << "Error : Trying to perform a MC Metropolis simulation without any move selection !" << endl;
        exit(-30);
    }

    // for keeping trace of move trials and acceptance for each movetype
    vector<int> nmvTrial(nmvtyp, 0);
    vector<int> nmvAcc(nmvtyp, 0);

    int natom = ens.getN();
    int imvtyp = 0;
    int imvatm = 0;

    double etot = ff.getE();
    double eold = 0.0, enew = 0.0;
    double de = 0.0;

    const vector<int>& nMoveAt = mvlist.getNMoveAtm();
    const vector<MOVETYPE>& movetypList = mvlist.getMoveTypeList();
    const vector<int**>& moveAtomList = mvlist.getMoveAtomList();
//     const vector<int**>& moveBondList = mvlist.getMoveBondList();
    const vector<int**>& movePivotList = mvlist.getMovePivotList();

    double r[3] = {0.0, 0.0, 0.0};
    double rang = 0.0;
    //    cout << "nmvtyp nmvat : \t" << nmvtyp << '\t' << nMoveAt << endl;

    // for storing 
    vector < tuple<double, double, double >> crdbackup(natom, tuple<double, double, double>(0.0, 0.0, 0.0));

    // MC metropolis main loop
    for ( int st = 1; st <= nsteps; st++ )
    {
        isAccepted = false;

        // get a random movetype and moveinstance from this movetype
        imvtyp = rndIntCandidate(nmvtyp); // get an int between 0 and (nmvtyp-1)
        imvatm = rndIntCandidate(nMoveAt[imvtyp]); // get an int between 0 and (nMoveAt[imvtyp]-1)

        // keep trace of number of trials for each mvtype
        nmvTrial[imvtyp]++;

        cout << "imvtyp : " << imvtyp << '\t' << "imvatm : " << imvatm << endl;

//         eold = E_moving_set(natom, nmvtyp, imvtyp, imvatm);
        eold = ff.getE();
        
        Atom::crd_backup_save(crdbackup, at_List, moveAtomList[imvtyp][imvatm]);

        // apply move
        switch ( movetypList[imvtyp] )
        {
            case TRN:
            {
                rndSphere(r);
                scaleVec(r, dmax);
                Move_TRN::translate_set(at_List, pbc, moveAtomList[imvtyp][imvatm],
                                        r[0], r[1], r[2]);
                break;
            }
            case ROT:
            {
                rndSphere(r);
                rang = dmax * rndUnifMove();
                Move_ROT::rotate_set(at_List, pbc, moveAtomList[imvtyp][imvatm],
                                     movePivotList[imvtyp][imvatm][0], rang, r);
                break;
            }
            default:
                break;
        }

//         enew = E_moving_set(natom, nmvtyp, imvtyp, imvatm);
        enew = ff.getE();
        
        de = enew - eold;
        apply_criterion(de);

        if ( isAccepted )
        {
            nmvAcc[imvtyp]++;
            etot += de;
        }
        else
        {
            Atom::crd_backup_load(crdbackup, at_List, moveAtomList[imvtyp][imvatm]);
        }

        //if necessary write trajectory
        if ( st % svFreq == 0 )
        {
            write_traj(st);
            //             fprintf(efile,"%d\t%lf\t%lf\n",st,etot,ff.getEtot());
            fprintf(efile, "%d\t%lf\n", st, etot);
        }

    } // end of MC metropolis main loop

    for ( int i = 0; i < nmvtyp; i++ )
    {
        cout << "For move type " << i << " : \n" << "TRIALS \t" << nmvTrial[i] << "\tACCEPTED \t" << nmvAcc[i];
        cout << "\tACCEPTANCE \t" << 100.0 * (double) nmvAcc[i] / (double) nmvTrial[i] << endl << endl;
    }

} // end of function run()

void MC_metropolis::apply_criterion(double de)
{
    const double dbl_epsilon = numeric_limits<double>::epsilon();

    cout << de << '\t';

    if ( de <= dbl_epsilon )
    {
        isAccepted = true;
    }
    else
    {
        double alpha = rndUnifAlpha();
        double beta = 1.0 / ((FField::rboltzui / FField::kcaltoiu) * ens.getTemp());
        double accf = exp(-beta * de);

        cout << alpha << '\t' << accf;

        if ( alpha < accf )
        {
            isAccepted = true;
        }
        else
        {
            isAccepted = false;
        }
    }

    cout << endl;
}
