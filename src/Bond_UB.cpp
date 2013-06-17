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

#include <iomanip>

#include "Bond_UB.h"

Bond_UB::Bond_UB() : Bond() { }

Bond_UB::Bond_UB(int _a1, int _a2, int _typ, double _k, double _r)
: Bond(_a1, _a2, _typ, _k, _r, 0.0) { }

Bond_UB::~Bond_UB() { }

std::ostream& operator<<(std::ostream& overloadStream, const Bond_UB& bnd_ub)
{
    bnd_ub.toString(overloadStream);

    return overloadStream;
}

void Bond_UB::toString(std::ostream& stream) const
{
    stream << std::fixed << std::setprecision(6);
    stream << "Urey-Bradley" << '\t';
    stream << at1 << '\t' << at2 << '\t' << type << '\t';
    stream << k << '\t' << r0 << '\t' << beta;
}

