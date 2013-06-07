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
#include <cstring>

#include <iostream>

#include "IO_CHARMM.h"

IO_CHARMM::IO_CHARMM(std::string configf_name ,std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens) : IO(_at_List,_pbc,_ens)
{
    conff = NULL;
    conff = fopen(configf_name.c_str(),"rt");
    if (conff==NULL)
    {
        std::cout << "Error while opening the coordinates files " << configf_name << std::endl;
        exit(-5);
    }
    
    read_CONF();
}

IO_CHARMM::~IO_CHARMM() {
}

void IO_CHARMM::read_CONF()
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
    
    if(lnatom != ens.getN())
    {
        std::cout << "Error : number of atoms at the top of coordinates file differs"
                     "from one from the input XML file."<< std::endl;
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

    fclose(conff);
}