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
#include <ctime>
#include <iostream>

#include "Atom.h"
#include "PerConditions.h"
#include "Ensemble.h"
#include "Ens_NVT.h"
#include "Ens_NPT.h"
#include "FField.h"
#include "MC.h"

int main()
{
    srand((unsigned int)time(0));

    // 1 : create vector of atoms
    int natom = 500;
    double T = 120.0;
    double boxL = 30.0;
    double boxAng = 90.0;
    int nsteps = 1000000;
    double dmax = 0.25;
    
    
    std::vector<Atom> lst;
    for ( int i = 0 ; i < natom ; i++ )
        lst.push_back( Atom(i,0) );

    // 2 : define periodic boundaries conditions
    pbcond pb = CUBIC;
    PerConditions pbc = PerConditions(pb);
    pbc.set_pbc_vectors(boxL);
    pbc.set_pbc_angles(boxAng);

    // 3 : define the Ensemble
    Ensemble ens = Ens_NVT(natom,pbc.computeVol(),T);
//    Ensemble ens = Ens_NPT(natom,300,T);

    // 4: ForceField
    FField ff = FField(lst,pbc,ens);

    // 5 : create and run MC simulation
    MC simulation = MC(lst,pbc,ens,ff,nsteps,dmax);
    simulation.run();
    
    return EXIT_SUCCESS;
}
