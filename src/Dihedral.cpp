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

#include <iomanip>

#include "Dihedral.h"

Dihedral::Dihedral()
{
}

Dihedral::Dihedral(int _a1, int _a2, int _a3, int _a4,
                   int _typ, int _ord, double _k, double _phi0, double _mult)
                   : at1(_a1), at2(_a2), at3(_a3), at4(_a4),
                   type(_typ), order(_ord), k(_k), phi0(_phi0), mult(_mult)
{
}

Dihedral::~Dihedral()
{
}

std::ostream& operator<<(std::ostream& overloadStream, const Dihedral& dihe)
{   
    dihe.toString(overloadStream);
    
    return overloadStream;
}

void Dihedral::toString(std::ostream& stream) const
{
    stream << std::fixed << std::setprecision(6);
    stream << "Dihedral" << '\t';
    stream << at1 << '\t' << at2 << '\t' << at3 << '\t'<< at4 << '\t' << type << '\t';
    stream << order << '\t' << k << '\t' << phi0 << '\t' << mult;
}

