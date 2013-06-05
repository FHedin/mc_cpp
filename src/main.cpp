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
#include <cstring>

#include <iostream>

#include "Parser.h"

#include "Atom.h"

#include "PerConditions.h"

#include "Ensemble.h"
#include "Ens_NVT.h"
#include "Ens_NPT.h"

#include "FField.h"

#include "MC.h"
#include "MC_metropolis.h"
#include "MC_spav.h"

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        cout << "Error with arguments processing : please provide the input file name : " << endl;
        cout << "Example : " << endl << argv[0] << " -i an_input_file.xml " << endl;
        exit(-1);            
    }
    
    //arguments parsing
    char *inpname = NULL;
    for (int i=1; i<argc; i++)
    {
        if (!strcmp(argv[i],"-i"))
            inpname = argv[++i];
        else
        {
            cout << "Error : Argument '" << argv[i] << "' not recognised. " << endl;
            exit(-2);
        }
    }

    // efficient xml parsing of parameters
    Parser_XML xmlfp(inpname);
    
    int natom =   xmlfp.val_from_attr<int>("N");
    double T =    xmlfp.val_from_attr<double>("T");
    
    double boxL =       xmlfp.val_from_attr<double>("a");
    double boxAng =     xmlfp.val_from_attr<double>("alpha");
    
    int nsteps =  xmlfp.val_from_attr<int>("nsteps");
    double dmax = xmlfp.val_from_attr<double>("dmax");
    
    int update_frequency = 100;
    
//    double we = 0.05;
//    int me = 1;
//    int ne = 25;
    
    std::vector<Atom> lst;
    for ( int i = 0 ; i < natom ; i++ )
        lst.push_back( Atom(i,0) );

    // 2 : define periodic boundaries conditions
    pbcond pb = CUBIC;
//    pbcond pb = NONE;
    PerConditions* pbc = new PerConditions(pb);
    pbc->set_pbc_vectors(boxL);
    pbc->set_pbc_angles(boxAng);
    
    // 3 : define the Ensemble
    Ensemble* ens = new Ens_NVT(natom,pbc->computeVol(),T);

    // 4: ForceField
    FField* ff = new FField(lst,*pbc,*ens);

    // 5 : create and run MC simulation
    MC* simulation = new MC_metropolis(lst,*pbc,*ens,*ff,nsteps,dmax,update_frequency);
//    MC* simulation = new MC_spav(lst,*pbc,*ens,*ff,nsteps,dmax,update_frequency,we,me,ne);
    simulation->run();
    
    /* freeing memory previously allocated with new */
    delete pbc;
    delete ens;
    delete ff;
    delete simulation;
    
    return EXIT_SUCCESS;
}

