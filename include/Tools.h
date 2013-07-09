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

#ifndef TOOLS_H
#define	TOOLS_H

#include <vector>

#include "Global_include.h"

#include "Atom.h"
#include "PerConditions.h"

class Tools
{
public:
    Tools();
    virtual ~Tools() = 0;

    static void str_rm_blank_spaces(std::string& str);

    static void str_to_lower_case(std::string& str);
    static void str_to_lower_case(char* str);

    static void vec_substract(const double a[3], const double b[3], double c[3]);

    static void getCentreOfMass(std::vector<Atom>& at_List, double cmass[3]);
    static void getCentreOfMass(std::vector<Atom>& at_List, int first, int last, double cmass[3]);
    static void getCentreOfMass(std::vector<Atom>& at_List, int moveAtmList[], double cmass[3]);

    static double distance2(const double a1[3], const double a2[3],
                            const PerConditions& pbc);
    static double distance2(const double a1[3], const double a2[3],
                            const PerConditions& pbc, double delta[3]);

    template <typename T>
    static T X2(T x);

    template <typename T>
    static T X3(T x);

    template <typename T>
    static T X4(T x);

    template <typename T>
    static T X6(T x);

    template <typename T>
    static T X12(T x);

private:

};

template <typename T>
T Tools::X2(T x)
{
    return (( x )*( x ) );
}

template <typename T>
T Tools::X3(T x)
{
    return (( x )*( x )*( x ) );
}

template <typename T>
T Tools::X4(T x)
{
    return (( x )*( x )*( x )*( x ) );
}

template <typename T>
T Tools::X6(T x)
{
    return (( x )*( x )*( x )*( x )*( x )*( x ) );
}

template <typename T>
T Tools::X12(T x)
{
    return (( x )*( x )*( x )*( x )*( x )*( x )*( x )*( x )*( x )*( x )*( x )*( x ) );
}


#endif	/* TOOLS_H */

