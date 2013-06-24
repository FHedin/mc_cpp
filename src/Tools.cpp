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

#include <cstring>

#include <string>
#include <algorithm> //for std::remove_if
#include <functional> //for std::ptr_fun

#include "Tools.h"

Tools::Tools() { }

Tools::~Tools() { }

//removes white spaces from input string

void Tools::str_rm_blank_spaces(std::string& str)
{

    str.erase(std::remove_if(str.begin(), str.end(), std::ptr_fun(isspace)), str.end());
}

// for comparisons, be sure that the string str is in low caps

void Tools::str_to_lower_case(std::string& str)
{
    std::transform(str.begin(), str.end(), str.begin(), std::ptr_fun(tolower));
}

void Tools::str_to_lower_case(char* str)
{
    size_t n = strlen(str);

    for ( size_t i = 0; i < n; i++ )
    {
        str[i] = (char) tolower(str[i]);
    }
}

// stores b-a in c

void Tools::vec_substract(const double a[3], const double b[3], double c[3])
{
    c[0] = b[0] - a[0];
    c[1] = b[1] - a[1];
    c[2] = b[2] - a[2];
}

void Tools::getCentreOfMass(std::vector<Atom>& at_List, double cmass[3])
{
    double crd[3];
    double mass, mtot = 0.0;
    cmass[0] = cmass[1] = cmass[2] = 0.0;

    for ( auto it : at_List )
    {
        it.getCoords(crd);
        mass = it.getMass();
        mtot += mass;
        cmass[0] += mass * crd[0];
        cmass[1] += mass * crd[1];
        cmass[2] += mass * crd[2];
    }
    cmass[0] /= mtot;
    cmass[1] /= mtot;
    cmass[2] /= mtot;
}

void Tools::getCentreOfMass(std::vector<Atom>& at_List, int first, int last, double cmass[3])
{
    double crd[3];
    double mass, mtot = 0.0;
    cmass[0] = cmass[1] = cmass[2] = 0.0;

    for ( int i = first; i <= last; i++ )
    {
        at_List[i].getCoords(crd);
        mass = at_List[i].getMass();
        mtot += mass;
        cmass[0] += mass * crd[0];
        cmass[1] += mass * crd[1];
        cmass[2] += mass * crd[2];
    }
    cmass[0] /= mtot;
    cmass[1] /= mtot;
    cmass[2] /= mtot;
}

void Tools::getCentreOfMass(std::vector<Atom>& at_List, int moveAtmList[], double cmass[3])
{
    double crd[3];
    double mass, mtot = 0.0;
    cmass[0] = cmass[1] = cmass[2] = 0.0;

    int ng = moveAtmList[0];
    int endng = ng + 2;
    int nn;
    int iaf, ial;
    for ( int it1 = 1; it1 <= ng; it1++ )
    {
        nn = moveAtmList[it1];
        for ( int it2 = endng; it2 <= nn; it2 += 2 )
        {
            iaf = moveAtmList[it2 - 1];
            ial = moveAtmList[it2];
            for ( int it3 = iaf; it3 <= ial; it3++ )
            {
                at_List[it3].getCoords(crd);
                mass = at_List[it3].getMass();
                mtot += mass;
                cmass[0] += mass * crd[0];
                cmass[1] += mass * crd[1];
                cmass[2] += mass * crd[2];
            }
        }
        endng = nn + 2;
    }

    cmass[0] /= mtot;
    cmass[1] /= mtot;
    cmass[2] /= mtot;
}

double Tools::distance2(const double a1[3], const double a2[3], const PerConditions& pbc)
{
    double r2;
    double delta[3];

    Tools::vec_substract(a1, a2, delta);

    pbc.applyPBC(delta);

    r2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];

    return r2;
}

double Tools::distance2(const double a1[3], const double a2[3], const PerConditions& pbc, double delta[3])
{
    double r2;

    Tools::vec_substract(a1, a2, delta);

    pbc.applyPBC(delta);

    r2 = delta[0] * delta[0] + delta[1] * delta[1] + delta[2] * delta[2];

    return r2;
}
