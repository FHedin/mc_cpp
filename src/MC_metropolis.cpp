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

#include "MC_metropolis.h"
#include "Move_TRN.h"
#include "Move_ROT.h"

using namespace std;

MC_metropolis::MC_metropolis(std::vector<Atom>& _at_List, PerConditions& _pbc,
                             Ensemble& _ens, FField& _ff, List_Moves& _mvlist, int _steps, double _dmax,
                             int _update_frequency) : MC(_at_List, _pbc, _ens, _ff, _mvlist)
{
    nsteps = _steps;
    dmax = _dmax;
    upFreq = _update_frequency;

    //assign random coordinates to the vector of atoms
    //Init();

    xyz = nullptr;
    xyz = fopen("tr.xyz", "w");
    //    efile = fopen("ener.dat", "w");
    //    pfile=fopen("press.dat","w");

    cout << "Initialising MC Metropolis simulation : found " << ens.getN() << " atoms. The ensemble is " << ens.whoami() << std::endl;
}

MC_metropolis::~MC_metropolis()
{
    fclose(xyz);
    //    fclose(efile);
    //    fclose(pfile);
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

    double efirst = ff.getEtot();
    const vector<int>& nMoveAt = mvlist.getNMoveAtm();
    const vector<MOVETYPE>& mvtypList = mvlist.getMoveTypeList();
    const vector<int**>& moveAtomList = mvlist.getMoveAtomList();
    const vector<int**>& movePivotList = mvlist.getMovePivotList();

    cout << "First Energy \t NMVTYP" << endl;
    cout << efirst << '\t' << nmvtyp << endl << endl;

    for ( int st = 1; st <= nsteps; st++ )
    {
        isAccepted = false;

        // get a random movetype and moveinstance from this movetype
        imvtyp = rndIntCandidate(nmvtyp); // get an int between 0 and (nmvtyp-1)
        imvatm = rndIntCandidate(nMoveAt[imvtyp]); // get an int between 0 and (nMoveAt[imvtyp]-1)

        //        cout << imvtyp << '\t' << imvatm << endl;

        // apply move
        switch ( mvtypList[imvtyp] )
        {
            case TRN:
            {
                double r[3];
                rndSphere(r);
                scaleVec(r, dmax);
                Move_TRN::translate_set(at_List, moveAtomList[imvtyp][imvatm],
                                        r[0], r[1], r[2]);
                break;
            }
            case ROT:
            {
                double r[3];
                rndSphere(r);
                double rang = 90.0 * rndUnifMove();
                Move_ROT::rotate_set(at_List, moveAtomList[imvtyp][imvatm],
                                     movePivotList[imvtyp][imvatm][0], rang, r);
                break;
            }
            default:
                break;
        }

        cout << "After step " << st << "\t new energy is : \t" << ff.getEtot() << endl;
        write_traj();

    } // MC metropolis main loop

    //    Atom newAt = at_List.at(0);
    //    int n = ens.getN();
    //    int candidate;
    //    double crd[3];
    //    int acc = 0, acc2 = 0;

    //    double e_tail=n*ff.tail_energy();
    //    double p_tail=ff.tail_pressure(); 

    //    std::cout << "Total Tail energy correction : " << e_tail << std::endl;
    //    std::cout << "Tail pressure correction : " << p_tail << std::endl;

    // initial energy and virial
    //    ens.setE(ff.getLJ(false));

    //    //equilibration cycle
    //    std::cout << "Step 1 : Equilibration cycle of " << 0.1 * nsteps << " steps." << std::endl;
    //    for ( int st = 1; st <= 0.1 * nsteps; st++ )
    //    {
    //        isAccepted = false;
    //
    //        n = ens.getN();
    //        candidate = rndCandidate(n);
    //
    //        Atom& oldAt = at_List.at(candidate);
    //
    //        oldAt.getCoords(crd);
    //        newAt.setCoords(crd);
    //
    //        move(newAt);
    //
    //        apply_criterion(oldAt, newAt, candidate);
    //
    //        if ( isAccepted )
    //        {
    //            newAt.getCoords(crd);
    //            oldAt.setCoords(crd);
    //            acc++;
    //        }
    //
    //        if ( st % upFreq == 0 )
    //        {
    //            adj_dmax(acc, upFreq);
    //            acc = 0;
    //        }
    //    }
    //    std::cout << "Step 1 : End of equilibration cycle." << std::endl;
    //
    // initial coordinates
    //    recentre();
    //    write_traj();

    //    ens.setE(ff.getLJ(false));
    //    fprintf(efile, "%lf\n", ens.getE());

    //production cycle
    //    std::cout << "Step 2 : Production cycle of " << nsteps << " steps." << std::endl;
    //    for ( int st = 1; st <= nsteps; st++ )
    //    {
    //        isAccepted = false;
    //
    //        //        for(int i=0; i<ens.getN(); i++)
    //        //                std::cout << ff.getExtraE(i) << "\t" ;
    //        //        std::cout << std::endl << std::endl;
    //
    //        n = ens.getN();
    //        candidate = rndCandidate(n);
    //
    //        Atom& oldAt = at_List.at(candidate);
    //
    //        oldAt.getCoords(crd);
    //        newAt.setCoords(crd);
    //
    //        move(newAt);
    //
    //        apply_criterion(oldAt, newAt, candidate);

    //        std::vector<Atom> candidateVector(at_List);
    //        move(candidateVector);
    //        apply_criterion(candidateVector);

    //        if ( isAccepted )
    //        {
    //            //            at_List = candidateVector;
    //            newAt.getCoords(crd);
    //            oldAt.setCoords(crd);
    //            acc++;
    //        }
    //
    //        if ( st % upFreq == 0 )
    //        {
    //            recentre();
    //            write_traj();
    //
    //            acc2 += acc;
    //            adj_dmax(acc, upFreq);
    //            acc = 0;
    //
    //            //            double p = ff.PressFromVirial(upFreq);
    //
    //            //            fprintf(efile, "%lf\n", ens.getE()/*+e_tail*/);
    //            //            fprintf(pfile,"%lf\n",p+p_tail);
    //
    //        }
    //    }
    //    std::cout << "Step 2 : End of production cycle." << std::endl;
    //    std::cout << "Acceptance is (%): " << 100.0 * (double) acc2 / (double) nsteps << std::endl;
    //    std::cout << "Final dmax is :" << dmax << std::endl;
}

//void MC_metropolis::apply_criterion(Atom const& oldAt, Atom const& newAt, int candidate)
//{
//    //for storing energies
//    double de = 0., deextra = 0.;
//    double e1 = 0., e2 = 0.;
//    double extra1 = 0., extra2 = 0.;
//
//    // contributions of internal energy (deltaU) and work of volume change (Press*deltaV)
//    double u = 0., v = 0.;
//
//    double alpha, acc;
//    double beta = 1.0 / ((FField::rboltzui / FField::kcaltoiu) * ens.getTemp());
//
//    //    u = ff.getLJ(oldAt, candidate, false);
//    // TODO : PdeltaV
//    e1 = u + v;
//    //    extra1 = ff.getExtraE(candidate);
//
//    //    u = ff.getLJ(newAt, candidate, false);
//    // TODO : PdeltaV
//    e2 = u + v;
//    //    extra2 = ff.getExtraE(candidate);
//
//    de = e2 - e1;
//    //    deextra = extra2 - extra1;
//
//    alpha = rndUnifAlpha();
//
//    acc = exp(-beta * (de + deextra));
//
//    if ( alpha <= acc )
//    {
//        isAccepted = true;
//        //        ens.addE(de);
//    }
//
//    //    std::cerr << "de \t deextra \t alpha \t acc \t test : " << de << "\t" << deextra << "\t" << alpha << "\t" << acc << "\t" << (alpha<=acc) << std::endl;
//}
//
//void MC_metropolis::apply_criterion(std::vector<Atom>& candidateVector)
//{
//    //for storing energies
//    double de = 0., deextra = 0.;
//    double e1 = 0., e2 = 0.;
//    double extra1 = 0., extra2 = 0.;
//
//    // contributions of internal energy (deltaU) and work of volume change (Press*deltaV)
//    double u = 0., v = 0.;
//
//    double alpha, acc;
//    double beta = 1.0 / ((FField::rboltzui / FField::kcaltoiu) * ens.getTemp());
//
//    //    u = ens.getE();
//    // TODO : PdeltaV
//    e1 = u + v;
//    //    extra1 = ff.getExtraE();
//
//    //    u = ff.getLJ(candidateVector, false);
//    // TODO : PdeltaV
//    e2 = u + v;
//    //    extra2 = ff.getExtraE();
//
//    de = e2 - e1;
//    //    deextra = extra2 - extra1;
//
//    alpha = rndUnifAlpha();
//
//    acc = exp(-beta * (de + deextra));
//
//    if ( alpha <= acc )
//    {
//        isAccepted = true;
//        //        ens.setE(e2);
//    }
//}


