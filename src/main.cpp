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

#include "Tools.h"

#include "PerConditions.h"

#include "Ensemble.h"
#include "Ens_NVT.h"
#include "Ens_NPT.h"

#include "Atom.h"

#include "FField.h"

#include "List_Exclude.h"

#include "MC.h"
#include "MC_metropolis.h"
#include "MC_spav.h"

#include "IO_CHARMM.h"
#include "IO_MDBAS.h"
#include "FField_MDBAS.h"

void get_simul_params_from_file(Parser_XML* xmlfp, PerConditions** pbc, Ensemble** ens,
        std::vector<Atom>& lst, FField** ff, List_exclude** exlst, MC** simulation);

using namespace std;

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        cerr << "Error with arguments processing : please provide the input file name : " << endl;
        cerr << "Example : " << endl << argv[0] << " -i an_input_file.xml " << endl;
        exit(-1);
    }

    //cmd line arguments parsing
    char* inpname = nullptr;
    for (int i = 1; i < argc; i++)
    {
        if (!strcmp(argv[i], "-i"))
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
    std::vector<Atom> lst;
    FField* ff = nullptr;
    List_exclude* exlst = nullptr;
    MC* simulation = nullptr;

    // efficient xml parsing of parameters
    //    xmlfp = new Parser_XML(inpname,true);
    xmlfp = new Parser_XML(inpname);
    get_simul_params_from_file(xmlfp, &pbc, &ens, lst, &ff, &exlst, &simulation);

    // run simulation immediately as everything was parsed before
    // simulation->run();
   
//    cout << *ff;
    ff->getEtot();
    
    /* freeing memory previously allocated with new */
    delete xmlfp;
    delete pbc;
    delete ens;
    delete ff;
    delete exlst;
    delete simulation;

    return EXIT_SUCCESS;
}

void get_simul_params_from_file(Parser_XML* xmlfp, PerConditions** pbc, Ensemble** ens,
        std::vector<Atom>& lst, FField** ff, List_exclude** exlst, MC** simulation)
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
    if (!enstype.compare("nvt"))
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

    *ff = new FField_MDBAS(lst, **pbc, **ens);

    // Atom list + coordinates reading
    string atMode = xmlfp->val_from_attr<string>("at_list");
    Tools::str_rm_blank_spaces(atMode);
    Tools::str_to_lower_case(atMode);
    if (!atMode.compare("repeat"))
    {
        string symb = xmlfp->val_from_attr<string>("symbol");
        double q = xmlfp->val_from_attr<double>("charge");
        double lj_eps = xmlfp->val_from_attr<double>("lj_epsilon");
        double lj_sig = xmlfp->val_from_attr<double>("lj_sigma");
        for (int i = 0; i < natom; i++)
        {
            lst.push_back(Atom(i, symb));
            lst.at(i).setCharge(q);
            lst.at(i).setEpsilon(lj_eps);
            lst.at(i).setSigma(lj_sig);
        }
    }
    else if (!atMode.compare("file"))
    {
        string corname = xmlfp->val_from_attr<string>("at_file");
        if (is_mdbas)
        {
            io = new IO_MDBAS(corname, fffile, **ff, lst, **pbc, **ens);
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
    *exlst = new List_exclude(**ff, **ens);
    (*ff)->setExcl(**exlst);

//    int nsteps = xmlfp->val_from_attr<int>("nsteps");
//    double dmax = xmlfp->val_from_attr<double>("dmax_value");
//    int update_frequency = xmlfp->val_from_attr<int>("dmax_update");
//
//    *simulation = new MC_metropolis(lst, **pbc, **ens, **ff, nsteps, dmax, update_frequency);

    delete io;
}

