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

#include <cstdlib>
#include <cstring>

#include <iostream>

#include "Parser.h"

#include "Tools.h"

#include "PerConditions.h"

#include "Ensemble.h"
#include "Ens_NVT.h"
#include "Ens_NPT.h"

#include "Atom.h"

#include "FField.h"

#include "List_Exclude.h"
#include "List_Moves.h"

#include "MC.h"
#include "MC_metropolis.h"
#include "MC_spav.h"

#include "IO_CHARMM.h"
#include "IO_MDBAS.h"
#include "FField_MDBAS.h"


void get_simul_params_from_file(Parser_XML* xmlfp, PerConditions** pbc, Ensemble** ens,
                                std::vector<Atom>& atomList, FField** ff, List_Exclude** exlst,
                                List_Moves** mvlist, MC** simulation);

using namespace std;

int main(int argc, char* argv[])
{
    if ( argc < 3 )
    {
        cerr << "Error with arguments processing : please provide the input file name : " << endl;
        cerr << "Example : " << endl << argv[0] << " -i an_input_file.xml " << endl;
        exit(-1);
    }

    //cmd line arguments parsing
    char* inpname = nullptr;
    for ( int i = 1; i < argc; i++ )
    {
        if ( !strcmp(argv[i], "-i") )
            inpname = argv[++i];
        else
        {
            cerr << "Error : Argument '" << argv[i] << "' not recognised. " << endl;
            exit(-2);
        }
    }

    Parser_XML* xmlfp = nullptr;
    PerConditions* pbc = nullptr;
    Ensemble* ens = nullptr;
    vector<Atom> atomList;
    FField* ff = nullptr;
    List_Exclude* exlst = nullptr;
    List_Moves* mvList = nullptr;
    MC* simulation = nullptr;

    // efficient xml parsing of parameters
    xmlfp = new Parser_XML(inpname, false);
    get_simul_params_from_file(xmlfp, &pbc, &ens, atomList, &ff, &exlst, &mvList, &simulation);

    delete xmlfp;

    //cout << *mvList;

    // run simulation immediately as everything was parsed before
     simulation->run();

    /* freeing memory previously allocated with new */
    delete simulation;
    delete mvList;
    delete exlst;
    delete ff;
    delete ens;
    delete pbc;

    return EXIT_SUCCESS;
}

void get_simul_params_from_file(Parser_XML* xmlfp, PerConditions** pbc, Ensemble** ens,
                                std::vector<Atom>& atomList, FField** ff, List_Exclude** exlst,
                                List_Moves** mvlist, MC** simulation)
{
    // box and periodic boundary conditions
    string pbtype = xmlfp->val_from_attr<string>("pbctype");
    double a = xmlfp->val_from_attr<double>("a");
    double b = xmlfp->val_from_attr<double>("b");
    double c = xmlfp->val_from_attr<double>("c");
    double alpha = xmlfp->val_from_attr<double>("alpha");
    double beta = xmlfp->val_from_attr<double>("beta");
    double gamma = xmlfp->val_from_attr<double>("gamma");
    *pbc = new PerConditions(pbtype, a, b, c, alpha, beta, gamma);

    // informations about the Ensemble
    string enstype = xmlfp->val_from_attr<string>("enstype");
    int natom = xmlfp->val_from_attr<int>("N");
    double T = xmlfp->val_from_attr<double>("T");
    Tools::str_rm_blank_spaces(enstype);
    Tools::str_to_lower_case(enstype);
    if ( !enstype.compare("nvt") )
    {
        *ens = new Ens_NVT(natom, (*pbc)->computeVol(), T);
    }
    else
    {
        cerr << "Error : the current version only supports the 'NVT' ensemble, "
                "with a possibly infinite volume (i.e. no pbc)." << std::endl;
        exit(-3);
    }

    // Forcefield parameters  + opening of coordinates files
    IO* io = nullptr;
    string ffmode = xmlfp->val_from_attr<string>("ff_mode");
    string fffile = xmlfp->val_from_attr<string>("ff_file");
    Tools::str_rm_blank_spaces(ffmode);
    Tools::str_to_lower_case(ffmode);
    Tools::str_rm_blank_spaces(fffile);

    bool is_mdbas = !ffmode.compare("mdbas");
#ifdef CHARMM_EXPERIMENTAL
    bool is_charmm = !ffmode.compare("charmm");
#endif 

    *ff = new FField_MDBAS(atomList, **pbc, **ens);

    // Atom list + coordinates reading
    string atMode = xmlfp->val_from_attr<string>("at_list");
    Tools::str_rm_blank_spaces(atMode);
    Tools::str_to_lower_case(atMode);
    if ( !atMode.compare("repeat") )
    {
        string symb = xmlfp->val_from_attr<string>("symbol");
        double q = xmlfp->val_from_attr<double>("charge");
        double lj_eps = xmlfp->val_from_attr<double>("lj_epsilon");
        double lj_sig = xmlfp->val_from_attr<double>("lj_sigma");
        for ( int i = 0; i < natom; i++ )
        {
            atomList.push_back(Atom(i, symb));
            atomList.at(i).setCharge(q);
            atomList.at(i).setEpsilon(lj_eps);
            atomList.at(i).setSigma(lj_sig);
        }
    }
    else if ( !atMode.compare("file") )
    {
        string corname = xmlfp->val_from_attr<string>("at_file");
        if ( is_mdbas )
        {
            io = new IO_MDBAS(corname, fffile, **ff, atomList, **pbc, **ens);
        }
        else
        {
            cerr << "Error : only MDBAS forcefield and coordinates (COR) supported ." << std::endl;
            exit(-7);
        }
    }
    else
    {
        cerr << "Error : the current version only supports "
                "the 'repeat' mode, or 'file' mode for the atomlist" << std::endl;
        exit(-4);
    }

    //build exclude list and link it to ff
    *exlst = new List_Exclude(**ff, **ens);
    (*ff)->setExcl(**exlst);

    // selection list and moves list
    // <move move_type="TRN" move_mode="ATOM" >
    string mvtyp = xmlfp->val_from_attr<string>("move_type");
    string mvmode = xmlfp->val_from_attr<string>("move_mode");
    *mvlist = new List_Moves(atomList, **ff, (*ens)->getN());
    string smode = xmlfp->val_from_attr<string>("sel_mode");
    string selec = xmlfp->val_from_node<string>("selection");
    (*mvlist)->addNewMoveType(mvtyp, mvmode, smode, selec);
    (*ff)->setMcmvlst(**mvlist);

    int nsteps = xmlfp->val_from_attr<int>("nsteps");
    double dmax = xmlfp->val_from_attr<double>("dmax_value");
    int save_frequency = xmlfp->val_from_attr<int>("save_each");

    *simulation = new MC_metropolis(atomList, **pbc, **ens, **ff, **mvlist, nsteps, dmax, save_frequency);

    delete io;
}

