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

#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>

#include "MC.h"

MC::MC(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens, FField& _ff, int _steps, double _dmax) : at_List(_at_List),pbc(_pbc),ens(_ens),ff(_ff)
{
    nsteps = _steps;
    dmax = _dmax;
    //assign random coordinates to the vector of atoms
    randInit();

    xyz=fopen("tr.xyz","w");
    efile=fopen("ener.dat","w");
    pfile=fopen("press.dat","w");

    std::cout << "Creating system : found " << ens.getN() << " atoms. The ensemble is " << ens.whoami() << std::endl;
}

MC::~MC()
{
    fclose(xyz);
    fclose(efile);
    fclose(pfile);
}

void MC::randInit()
{
    double crd[3];
    double pbv[3];

    switch (pbc.getType())
    {
        case NONE:
            pbv[0] = ens.getN();
            pbv[1] = ens.getN();
            pbv[2] = ens.getN();
            break;
        default:
            pbc.get_pbc_vectors(pbv);
            break;
    }

    for(std::vector<Atom>::iterator it = at_List.begin() ; it != at_List.end() ; ++it)
    {
        crd[0] = pbv[0] * ( -1. + 2.*(double)rand()/RAND_MAX );
        crd[1] = pbv[1] * ( -1. + 2.*(double)rand()/RAND_MAX );
        crd[2] = pbv[2] * ( -1. + 2.*(double)rand()/RAND_MAX );

        it->setCoords(crd);
        pbc.applyPBC(*it);
    }
}

void MC::run()
{
    Atom newAt;
    int n;
    int candidate;
    double crd[3];
    int acc=0,acc2=0;

    double e_tail=ff.tail_energy();
    double p_tail=ff.tail_pressure(); 

    //equilibration cycle
    std::cout << "Step 1 : Equilibration cycle of " << 0.1*nsteps << " steps." << std::endl;
    for (int st=1 ; st <= 0.1*nsteps ; st++)
    {
        isAccepted = false;

        n = ens.getN() ;
        candidate = (int) n*((double)rand()/RAND_MAX);

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
        
        if (st%1000 == 0)
        {
            adj_dmax(acc,1000);
            acc=0;
        }
    }
    std::cout << "Step 1 : End of equilibration cycle." << std::endl;

    // initial coordinates
    write_traj();
    // initial energy and virial
    ff.getLJV(false);
    ens.setE(ff.get_E());
    ff.resetE();
    ff.resetW();

    //production cycle
    std::cout << "Step 2 : Production cycle of " << nsteps << " steps." << std::endl;
    for (int st=1 ; st <= nsteps ; st++)
    {
        isAccepted = false;

        n = ens.getN() ;
        candidate = (int) n*((double)rand()/RAND_MAX);

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

        if (st%1000 == 0)
        {
            write_traj();
            
            acc2 += acc;
            adj_dmax(acc,1000);
            acc=0;
            
            double e = ens.getE() + e_tail ;
            double p = ff.PressFromVirial(1000) + p_tail ;
            
            fprintf(efile,"%lf\n",e);
            fprintf(pfile,"%lf\n",p);
            
            ff.resetW();
        }
    }
    std::cout << "Step 2 : End of production cycle." << std::endl;
    std::cout << "Acceptance is (%): " << 100.0*(double)acc2/(double)nsteps << std::endl;
    std::cout << "Final dmax is :" << dmax << std::endl;
}

void MC::move(Atom& newAt) const
{

    double initial[3]={0.} , trial[3]={0.} ;

    newAt.getCoords(initial);

    trial[0] = initial[0] + dmax*(-1.+2.*(double)rand()/RAND_MAX);
    trial[1] = initial[1] + dmax*(-1.+2.*(double)rand()/RAND_MAX);
    trial[2] = initial[2] + dmax*(-1.+2.*(double)rand()/RAND_MAX);

    newAt.setCoords(trial);

    pbc.applyPBC(newAt);
}

//random move for the whole box
void MC::move()
{

    double initial[3]={0.} , trial[3]={0.} ;

    for(std::vector<Atom>::iterator it = at_List.begin() ; it != at_List.end() ; ++it)
    {
        it->getCoords(initial);

        trial[0] = initial[0] + dmax*(-1.+2.*(double)rand()/RAND_MAX);
        trial[1] = initial[1] + dmax*(-1.+2.*(double)rand()/RAND_MAX);
        trial[2] = initial[2] + dmax*(-1.+2.*(double)rand()/RAND_MAX);

        it->setCoords(trial);

        pbc.applyPBC(*it);
    }
}

void MC::apply_criterion(Atom const& oldAt, Atom const& newAt, int candidate)
{
    //for storing energies
    double de=0. , e1=0. , e2=0. ;
    // contributions of internal energy (deltaU) and work of volume change (Press*deltaV)
    double u=0. , v=0. ;

    double alpha , acc ;

    double beta = 1.0/(FField::kb_ch*ens.getTemp()) ;

    ff.getLJV(oldAt,candidate,false);
    u = ff.get_E();
    // TODO : PdeltaV
    e1 = u + v;
    ff.resetE();

    ff.getLJV(newAt,candidate,false);
    u = ff.get_E();
    // TODO : PdeltaV
    e2 = u + v;
    ff.resetE();

    de = e2 - e1;

    alpha = (double)rand()/RAND_MAX;

    acc = exp(-beta*(de));

    if ( alpha <= acc)
    {
        isAccepted = true;
        ens.addE(de);
    }
}

void MC::write_traj() const
{
    double crd[3] = {0.};
    int n = ens.getN();

    fprintf(xyz,"%d\n",n);
    fprintf(xyz,"\n");

    for (std::vector<Atom>::iterator it = at_List.begin() ; it != at_List.end() ; ++it)
    {
        it->getCoords(crd);
        fprintf(xyz,"Ar\t%10.5lf\t%10.5lf\t%10.5lf\n",crd[0],crd[1],crd[2]);
    }
}

void MC::adj_dmax(double acc, double each)
{
//    std::cout << "dmax update : " << dmax << " --> "; 
    (acc/each)<=0.5 ? dmax*=0.95 : dmax*=1.05;
//    std::cout << dmax << " : targeting acceptance of 50 % " << std::endl;
    
    double pbv[3];
    pbc.get_pbc_vectors(pbv);
    
//    dmax > (pbv[0]/2.0 - 1.0) ? dmax = pbv[0]/2.0 - 1.0  : dmax;
//    dmax < 0.5 ? dmax = 0.5 : dmax;
}

