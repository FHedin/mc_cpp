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

#include <iostream>

#include "Atom.h"

#include "PerConditions.h"

#include "Ensemble.h"
#include "Ens_NVT.h"
#include "Ens_NPT.h"

#include "FField.h"

#include "MC.h"
#include "MC_metropolis.h"
#include "MC_spav.h"

int main()
{
    // 1 : create vector of atoms
    int natom = 55;
    double T = 0.25;
//    double boxL = 10.0;
//    double boxAng = 90.0;
    int nsteps = 100000;
    double dmax = 0.10;
    int update_frequency = 10000;
    
    double we = 0.05;
    int me = 1;
    int ne = 25;
    
    std::vector<Atom> lst;
    for ( int i = 0 ; i < natom ; i++ )
        lst.push_back( Atom(i,0) );

    // 2 : define periodic boundaries conditions
//    pbcond pb = CUBIC;
//    pbc.set_pbc_vectors(boxL);
//    pbc.set_pbc_angles(boxAng);
    
    pbcond pb = NONE;
    PerConditions* pbc = new PerConditions(pb);

    // 3 : define the Ensemble
    Ensemble* ens = new Ens_NVT(natom,pbc->computeVol(),T);

    // 4: ForceField
    FField* ff = new FField(lst,*pbc,*ens);

    // 5 : create and run MC simulation
//    MC* simulation = new MC_metropolis(lst,*pbc,*ens,*ff,nsteps,dmax,update_frequency);
    MC* simulation = new MC_spav(lst,*pbc,*ens,*ff,nsteps,dmax,update_frequency,we,me,ne);
    simulation->run();
    
    /* freeing memory previously allocated with new */
    delete pbc;
    delete ens;
    delete ff;
    delete simulation;
    
    return 0;
}

