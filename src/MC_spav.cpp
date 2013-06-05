/*
 *  mc_cpp : A basic Monte Carlo simulations software.
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

#ifdef SPAV_EXPERIMENTAL

#include <cstdio>
#include <iostream>

#include "MC_spav.h"

MC_spav::MC_spav(std::vector<Atom>& _at_List, PerConditions& _pbc, 
        Ensemble& _ens, FField& _ff, int _steps, double _dmax, 
        int _update_frequency, double _we, int _me, int _ne) : MC(_at_List,_pbc,_ens,_ff)
{
    nsteps = _steps;
    dmax = _dmax;
    upFreq = _update_frequency;
    
    we = _we;
    me = _me;
    ne = _ne;
    
    distributionNormal = std::normal_distribution<double>(0.0,we);
    
    oldAtList.resize(ne);
    newAtList.resize(ne);
    oldEList.resize(ne,0.0);
    newEList.resize(ne,0.0);
//    Sold.resize(me,0.0);
//    Snew.resize(me,0.0);
    deltaM.resize(me,0.0);
    
    //assign random coordinates to the vector of atoms
    Init();
    
    xyz=fopen("tr.xyz","w");
    efile=fopen("ener.dat","w");
//    pfile=fopen("press.dat","w");
    
    std::cout << "Initialising MC spatial averaging simulation : found " << ens.getN() << " atoms. The ensemble is " << ens.whoami() << std::endl;
    std::cout << "Spatial averaging parameters : we = " << we << "\t ne = " << ne << std::endl;
}

MC_spav::~MC_spav()
{
    fclose(xyz);
    fclose(efile);
//    fclose(pfile);
}

void MC_spav::run()
{
    Atom newAt;
    int n;
    int candidate;
    double crd[3];
    int acc=0,acc2=0;
    
    //equilibration cycle
    /*
     *  TODO
     */
    
    // initial coordinates
    write_traj();
    // initial energy and virial
    ens.setE(ff.getLJV(false));
    
    fprintf(efile,"%lf\n",ens.getE());
    
    //production cycle
    std::cout << "Step 2 : Production cycle of " << nsteps << " steps." << std::endl;
    for (int st=1 ; st <= nsteps ; st++)
    {
        isAccepted = false;
        
        n = ens.getN() ;
        candidate = rndCandidate(n);
        
//        std::cerr << "candidate is " << candidate << std::endl;
        
        Atom& oldAt = at_List.at(candidate);
        oldAt.getCoords(crd);
        newAt.setCoords(crd);
        move(newAt);
        
        apply_criterion(oldAt,newAt,candidate);
        
        if(isAccepted)
        {
            newAt.getCoords(crd);
            oldAt.setCoords(crd);
            acc++;
        }
        
        if (st%upFreq == 0)
        {
            recentre();
            write_traj();
            
            acc2 += acc;
            adj_dmax(acc,upFreq);
            acc=0;
            
            double e = ens.getE() ;
            
            fprintf(efile,"%.15lf\n",e);
        }
    }
    std::cout << "Step 2 : End of production cycle." << std::endl;
    std::cout << "Acceptance is (%): " << 100.0*(double)acc2/(double)nsteps << std::endl;
    std::cout << "Final dmax is :" << dmax << std::endl;
    
}

void MC_spav::buildSpavConfigs(Atom const& oldAt, Atom const& newAt)
{
    double add[3];
    
//    for(Atom& ref : oldAtList)
//        ref = oldAt;
//    
//    for(Atom& ref : newAtList)
//        ref = newAt;
    
//    std::cerr << oldAt.getX() << " " << oldAt.getY() << " " << oldAt.getZ() << std::endl;
//    std::cerr << newAt.getX() << " " << newAt.getY() << " " << newAt.getZ() << std::endl;
//    std::cerr << std::endl;
    
    oldAtList.assign(ne,oldAt);
    newAtList.assign(ne,newAt);
    
    for(Atom& ref : oldAtList)
    {
        rndNorm(add);
        ref.addCoords(add);
//        std::cerr << ref.getX() << " " << ref.getY() << " " << ref.getZ() << std::endl;
    }
//    std::cerr << std::endl;
    
    for(Atom& ref : newAtList)
    {
        rndNorm(add);
        ref.addCoords(add);
//        std::cerr << ref.getX() << " " << ref.getY() << " " << ref.getZ() << std::endl;
    }
//    std::cerr << std::endl << std::endl << std::endl;
}

void MC_spav::apply_criterion(Atom const& oldAt, Atom const& newAt, int candidate)
{
    //for storing energies
    double de=0., deextra=0. ;
    double e1=0., e2=0. ;
    double extra1=0., extra2=0.;
    
    // contributions of internal energy (deltaU) and work of volume change (Press*deltaV)
    double u=0. , v=0. ;
    
    // other variables
    double alpha , mcacc , spacc;
    double beta = 1.0/(FField::kb_ch*ens.getTemp()) ;
        
    /* ------------------------------------------------------------------- */
    // Metropolis like
    u = ff.getLJV(oldAt,candidate,false);
    // TODO : PdeltaV
    e1 = u + v;
//    extra1 = ff.getExtraE(candidate);

    u = ff.getLJV(newAt,candidate,false);
    // TODO : PdeltaV
    e2 = u + v;
//    extra2 = ff.getExtraE(candidate);

    de = e2 - e1;
//    deextra = extra2 - extra1;
    
    alpha = rndUnifAlpha();
    
    mcacc = exp(-beta*(de+deextra));
    
    std::cerr << "MC : DE \t alpha \t mcacc : " << de+deextra << "\t" << alpha << "\t" << mcacc << std::endl;
    
    if ( alpha <= mcacc)
    {
        isAccepted = true;
        ens.addE(de);
    }
    else
    {

        /* ------------------------------------------------------------------- */
        // spatial averaging 
        double Sold=0.0 , Snew=0.0;
        double delta = 0.0;

        for(int i=0; i<me; i++)
        {
            buildSpavConfigs(oldAt,newAt);

            for(int j=0; j<ne; j++)
            {
                oldEList.at(j) = ff.getLJV(oldAtList.at(j),candidate,false);
//                oldEList.at(j)+= ff.getExtraE(candidate);

                newEList.at(j) = ff.getLJV(newAtList.at(j),candidate,false);
//                newEList.at(j)+= ff.getExtraE(candidate);
            }

    //    //        std::cerr << "Old BolzE : " << std::endl;
    //        for(double& ref : oldEList)
    //        {
    //            ref = exp(-beta*ref);
    //    //            std::cerr << ref << "\t";
    //        }
    //    //        std::cerr << std::endl;
    //
    //    //        std::cerr << "New BolzE : " << std::endl;
    //        for(double& ref : newEList)
    //        {
    //            ref = exp(-beta*ref);
    //    //            std::cerr << ref << "\t";
    //        }
    //    //        std::cerr << std::endl;

//            Sold.assign(ne,0.0);
//            Snew.assign(ne,0.0);

            for(double& ref : oldEList)
                Sold += ref;

            for(double& ref : newEList)
                Snew += ref;

        //    for(int i=0; i<ne; i++)
        //    {
        //        Sold += oldEList.at(i);
        //        Snew += newEList.at(i);
        //    }

            deltaM.at(i) = -1.0*log(Snew/Sold);
        }
        
        //        spacc = exp(-beta*spacc);

            std::cerr << "SP : alpha \t spacc : " << alpha << "\t" << spacc << std::endl;

            if ( alpha <= spacc )
            {
                isAccepted = true;
                ens.addE(de);
            }
//        }
        std::cerr << std::endl;
    }
}

void MC_spav::apply_criterion(std::vector<Atom>& candidateVector)
{
    
}

double MC_spav::rndNorm()
{
    return distributionNormal(generator);
}

void MC_spav::rndNorm(double _crd[3])
{
    _crd[0] = distributionNormal(generator);
    _crd[1] = distributionNormal(generator);
    _crd[2] = distributionNormal(generator);
}


#endif //SPAV_EXPERIMENTAL