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

#include "Bond.hpp"

Bond::Bond()
{
    at1 = at2 = type = -1;
    k = r0 = beta = -1.0;
}

Bond::Bond(int _a1, int _a2, int _typ, double _k,
           double _r, double _beta) : at1(_a1), at2(_a2), type(_typ), k(_k), r0(_r), beta(_beta) { }

Bond::~Bond() { }

double Bond::getBeta() const
{
    return beta;
}

double Bond::getR0() const
{
    return r0;
}

double Bond::getK() const
{
    return k;
}

int Bond::getType() const
{
    return type;
}

int Bond::getAt2() const
{
    return at2;
}

int Bond::getAt1() const
{
    return at1;
}

std::ostream& operator<<(std::ostream& overloadStream, const Bond& bnd)
{
    bnd.toString(overloadStream);

    return overloadStream;
}

void Bond::toString(std::ostream& stream) const
{
    stream << std::fixed << std::setprecision(6);
    stream << "Bond" << '\t';
    stream << at1 << '\t' << at2 << '\t' << type << '\t';
    stream << k << '\t' << r0 << '\t' << beta;
}
