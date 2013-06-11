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

#include "Bond.h"

Bond::Bond()
{
    at1=at2=type=-1;
    k=r0=beta=-1.0;
}

Bond::Bond(int _a1, int _a2, int _typ, double _k,
        double _r, double _beta) : at1(_a1), at2(_a2), type(_typ), k(_k), r0(_r), beta(_beta)
{
}

Bond::~Bond()
{
}

std::ostream& operator<<(std::ostream& overloadStream, const Bond& bnd)
{
    overloadStream << std::fixed << std::setprecision(6);
    overloadStream << "Bond" << '\t';
    overloadStream << bnd.at1 << '\t' << bnd.at2 << '\t' << bnd.type << '\t';
    overloadStream << bnd.k << '\t' << bnd.r0 << '\t' << bnd.beta;
    
    return overloadStream;
}