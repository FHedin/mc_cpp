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

#include <cmath>
#include <cstring>

#include <iostream>
#include <vector>

#include "IO_MDBAS.h"

#include "Tools.h"

#include "Atom.h"

#include "FField.h"

#include "Bond.h"
#include "Bond_UB.h"

#include "Angle.h"

#include "Dihedral.h"
#include "Dihedral_improper.h"
#include "Dihedral_improper.h"

IO_MDBAS::IO_MDBAS(std::string configf_name, std::string forfieldf_name,
        std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens) : IO(_at_List, _pbc, _ens)
{
    conff = NULL;
    forff = NULL;

    conff = fopen(configf_name.c_str(), "rt");
    if (conff == NULL)
    {
        std::cerr << "Error while opening the coordinates files " << configf_name << std::endl;
        exit(-5);
    }

    read_coord();
    fclose(conff);

    forff = fopen(forfieldf_name.c_str(), "rt");
    if (conff == NULL)
    {
        std::cerr << "Error while opening the coordinates files " << configf_name << std::endl;
        exit(-8);
    }

    read_ff();
    fclose(forff);
}

IO_MDBAS::~IO_MDBAS()
{
}

void IO_MDBAS::read_coord()
{
    char buff1[1024] = "", *buff2 = NULL;

    char ren[5] = "", atl[5] = "", sen[5] = "";
    int atn, res, ire;
    double wei, xx, yy, zz;
    int lnatom;

    while (fgets(buff1, 1024, conff) != NULL)
    {
        if (buff1[0] != '*')
            break;
    }

    buff2 = strtok(buff1, " \n\t");
    lnatom = atoi(buff2);

    if (lnatom != ens.getN())
    {
        std::cerr << "Error : number of atoms at the top of coordinates file differs"
                "from one from the input XML file." << std::endl;
        exit(-6);
    }

    at_List.resize(lnatom);

    for (int i = 0; i < lnatom; i++)
    {
        fscanf(conff, "%d %d %4s %4s %lf %lf %lf %4s %d %lf", &atn, &ire, ren, atl, &xx, &yy, &zz, sen, &res, &wei);

        at_List.at(i).setId(atn);

        at_List.at(i).setX(xx);
        at_List.at(i).setY(yy);
        at_List.at(i).setZ(zz);

        at_List.at(i).setResidue_id_global(ire);
        at_List.at(i).setResidue_id_seg(res);

        at_List.at(i).setSymbol(atl);
        at_List.at(i).setRes_label(ren);
        at_List.at(i).setSeg_label(sen);
    }
}

void IO_MDBAS::read_ff()
{
    char buff1[1024] = "", *buff2 = NULL, buff3[1024] = "", *buff4 = NULL;
    int i, k, ia, ib;
    int nAtom;

    int nBond = 0;
    int nConst = 0;
    int nUb = 0;
    int nAngle = 0;
    int nDihedral = 0;
    int nImproper = 0;

    std::vector<Bond> bndList;
    std::vector<Bond_UB> ubList;
    std::vector<Angle> angList;
    std::vector<Dihedral> diheList;
    std::vector<Dihedral_improper> imprList;

    while (fgets(buff1, 1024, forff) != NULL)
    {
        Tools::str_to_lower_case(buff1);

        if (buff1[0] == '#')
            continue;

        buff2 = strtok(buff1, " \n\t");

        if (!strcmp(buff2, "atoms"))
        {

            nAtom = atoi(strtok(NULL, " \n\t"));
            if (nAtom != ens.getN())
            {
                std::cerr << "Error : number of atoms at the top of forcefield file differs"
                        "from one from the input XML file." << std::endl;
                exit(-9);
            }

            k = 0;
            while (k < nAtom)
            {
                if (fgets(buff3, 1024, forff) != NULL)
                {
                    if (buff3[0] == '#')
                        continue;

                    i = atoi(strtok(buff3, " \n\t")) - 1;
                    //                    if (i != k)
                    //                        my_error(PSF_FILE_ERROR, __FILE__, __LINE__, 0);

                    buff4 = strtok(NULL, " \n\t");

                    at_List.at(i).setType(atoi(strtok(NULL, " \n\t")));
                    at_List.at(i).setCharge(atof(strtok(NULL, " \n\t")));
                    at_List.at(i).setMass(atof(strtok(NULL, " \n\t")));
                    at_List.at(i).setIs_frozen((bool) atoi(strtok(NULL, " \n\t")));

                    k++;
                }
                else
                {
                    std::cerr << "Error when reading forcefield file : check the source file :"
                            << __FILE__ << "\t at line : " << __LINE__ << std::endl;
                    exit(-10);
                }
            }
        } // end of if (!strcmp(buff2, "atoms"))
        else if (!strcmp(buff2, "bonds"))
        {
            nBond = atoi(strtok(NULL, " \n\t"));

            if (nBond == 0)
                continue;

//            bndList.resize(nBond);

            int a, b, type;
            double kst, r0, beta;

            k = 0;
            while (k < nBond)
            {
                if (fgets(buff3, 1024, forff) != NULL)
                {
                    if (buff3[0] == '#')
                        continue;

                    a = atoi(strtok(buff3, " \n\t")) - 1;
                    b = atoi(strtok(NULL, " \n\t")) - 1;
                    type = atoi(strtok(NULL, " \n\t"));
                    kst = atof(strtok(NULL, " \n\t"))* FField::kcaltoiu;
                    r0 = atof(strtok(NULL, " \n\t"));
                    beta = atof(strtok(NULL, " \n\t"));

                    bndList.push_back(Bond(a, b, type, kst, r0, beta));

                    k++;
                }
                else
                {
                    std::cerr << "Error when reading forcefield file : check the source file :"
                            << __FILE__ << "\t at line : " << __LINE__ << std::endl;
                    exit(-10);
                }
            }
        } // end of else if (!strcmp(buff2, "bonds"))
        else if (!strcmp(buff2, "constraints"))
        {
            std::cerr << "Warning : constraints not implemented for the moment. Skipping section ... " << std::endl;

            nConst = atoi(strtok(NULL, " \n\t"));
            if (nConst == 0)
                continue;

            k = 0;
            while (k < nConst)
            {
                if (fgets(buff3, 1024, forff) != NULL)
                {
                    if (buff3[0] == '#')
                        continue;
                }
                else
                {

                }
            }

        } // end of else if (!strcmp(buff2, "constraints"))
        else if (!strcmp(buff2, "urey-bradley"))
        {
            nUb = atoi(strtok(NULL, " \n\t"));

            if (nUb == 0)
                continue;

//            ubList.resize(nUb);

            int a, b, type;
            double kst, r0;

            k = 0;
            while (k < nUb)
            {
                if (fgets(buff3, 1024, forff) != NULL)
                {
                    if (buff3[0] == '#')
                        continue;

                    a = atoi(strtok(buff3, " \n\t")) - 1;
                    b = atoi(strtok(NULL, " \n\t")) - 1;
                    type = atoi(strtok(NULL, " \n\t"));
                    kst = atof(strtok(NULL, " \n\t")) * FField::kcaltoiu;
                    r0 = atof(strtok(NULL, " \n\t"));

                    bndList.push_back(Bond_UB(a, b, type, kst, r0));

                    k++;
                }
                else
                {
                    std::cerr << "Error when reading forcefield file : check the source file :"
                            << __FILE__ << "\t at line : " << __LINE__ << std::endl;
                    exit(-10);
                }
            }
        } // end of else if (!strcmp(buff2, "urey-bradley"))
        else if (!strcmp(buff2, "angles"))
        {
            nAngle = atoi(strtok(NULL, " \n\t"));

            if (nAngle == 0)
                continue;

//            angList.resize(nAngle);

            int a, b, c, type;
            double kst, theta0;

            k = 0;
            while (k < nAngle)
            {
                if (fgets(buff3, 1024, forff) != NULL)
                {
                    if (buff3[0] == '#')
                        continue;

                    a = atoi(strtok(buff3, " \n\t")) - 1;
                    b = atoi(strtok(NULL, " \n\t")) - 1;
                    c = atoi(strtok(NULL, " \n\t")) - 1;

                    type = atoi(strtok(NULL, " \n\t"));

                    kst = atof(strtok(NULL, " \n\t")) * FField::kcaltoiu;
                    theta0 = atof(strtok(NULL, " \n\t")) * FField::PI / 180.;

                    angList.push_back(Angle(a, b, c, type, kst, theta0));

                    k++;
                }
                else
                {
                    std::cerr << "Error when reading forcefield file : check the source file :"
                            << __FILE__ << "\t at line : " << __LINE__ << std::endl;
                    exit(-10);
                }
            }
        } // end of else if (!strcmp(buff2, "angles"))
        else if (!strcmp(buff2, "dihedrals"))
        {
            nDihedral = atoi(strtok(NULL, " \n\t"));

            if (nDihedral == 0)
                continue;

//            diheList.resize(nDihedral);

            int a, b, c, d, type, order;
            double kst, phi0, mult;

            k = 0;
            while (k < nDihedral)
            {
                if (fgets(buff3, 1024, forff) != NULL)
                {
                    if (buff3[0] == '#')
                        continue;
                    
//                    std::cout << buff3;

                    a = atoi(strtok(buff3, " \n\t")) - 1;
                    b = atoi(strtok(NULL, " \n\t")) - 1;
                    c = atoi(strtok(NULL, " \n\t")) - 1;
                    d = atoi(strtok(NULL, " \n\t")) - 1;

                    type = atoi(strtok(NULL, " \n\t"));
                    order = atoi(strtok(NULL, " \n\t"));

                    kst = atof(strtok(NULL, " \n\t")) * FField::kcaltoiu;
                    phi0 = atof(strtok(NULL, " \n\t")) * FField::PI / 180.;
                    mult = atof(strtok(NULL, " \n\t"));

                    diheList.push_back(Dihedral(a, b, c, d, type, order, kst, phi0, mult));

                    k++;
                }
                else
                {
                    std::cerr << "Error when reading forcefield file : check the source file :"
                            << __FILE__ << "\t at line : " << __LINE__ << std::endl;
                    exit(-10);
                }
            }
        } // end of else if (!strcmp(buff2, "dihedrals"))
        else if (!strcmp(buff2, "impropers"))
        {
            nImproper = atoi(strtok(NULL, " \n\t"));

            if (nImproper == 0)
                continue;

//            imprList.resize(nImproper);

            int a, b, c, d, type, order;
            double kst, phi0, mult;

            k = 0;
            while (k < nImproper)
            {
                if (fgets(buff3, 1024, forff) != NULL)
                {
                    if (buff3[0] == '#')
                        continue;
                    
                    a = atoi(strtok(buff3, " \n\t")) - 1;
                    b = atoi(strtok(NULL, " \n\t")) - 1;
                    c = atoi(strtok(NULL, " \n\t")) - 1;
                    d = atoi(strtok(NULL, " \n\t")) - 1;

                    type = atoi(strtok(NULL, " \n\t"));
                    order = atoi(strtok(NULL, " \n\t"));

                    kst = atof(strtok(NULL, " \n\t")) * FField::kcaltoiu;
                    phi0 = atof(strtok(NULL, " \n\t")) * FField::PI / 180.;
                    mult = atof(strtok(NULL, " \n\t"));
                    
                    imprList.push_back(Dihedral_improper(a,b,c,d,type,order,kst,phi0,mult));
                    
                    k++;
                }
                else
                {
                    std::cerr << "Error when reading forcefield file : check the source file :"
                            << __FILE__ << "\t at line : " << __LINE__ << std::endl;
                    exit(-10);
                }
            }
        } // end of !strcmp(buff2, "impropers")
        else if (!strcmp(buff2, "vdw"))
        {
            nAtom = atoi(strtok(NULL, " \n\t"));
            if (nAtom != ens.getN())
            {
                std::cerr << "Error : number of vdw parameters is not the same "
                             "that the number of atoms of the system" << std::endl;
                exit(-11);
            }
            
            int type;
            double bet;
            
            k = 0;
            while (k < nAtom)
            {
                if (fgets(buff3, 1024, forff) != NULL)
                {
                    if (buff3[0] == '#')
                        continue;
                    
                    i = atoi(strtok(buff3, " \n\t")) - 1;
                    if (i != k)
                    {
                        std::cerr << "Error : Atom missing for vdw parameters :"
                        " atom " << i << " not available in forcefield file." << std::endl;
                        exit(-13);
                    }
                    
                    type = atoi(strtok(NULL, " \n\t"));
                    
                    at_List.at(i).setEpsilon( sqrt(atof(strtok(NULL, " \n\t")) * FField::kcaltoiu) );
                    at_List.at(i).setSigma( atof(strtok(NULL, " \n\t")) );
                    bet = atof(strtok(NULL, " \n\t"));
                    
                    at_List.at(i).setEpsilon14( sqrt(atof(strtok(NULL, " \n\t")) * FField::kcaltoiu) );
                    at_List.at(i).setSigma14( atof(strtok(NULL, " \n\t")) );
                    bet = atof(strtok(NULL, " \n\t"));
                    
                    k++;
                }
                else
                {
                    std::cerr << "Error when reading forcefield file : check the source file :"
                            << __FILE__ << "\t at line : " << __LINE__ << std::endl;
                    exit(-10); 
                }
            }
        } // end of !strcmp(buff2, "vdw")
        else if (!strcmp(buff2, "end"))
        {
            break;
        }
        else
        {
            std::cerr << "Error when reading forcefield file : unknown keyword '" << buff2 << "' : "
                      << __FILE__ << "\t at line : " << __LINE__ << std::endl;
            exit(-12); 
        }
    } //end of while (fgets(buff1, 1024, forff) != NULL)

    for(i=0 ; i<nBond; i++ )
        std::cout << bndList.at(i) << std::endl ;
    
} // end of function

