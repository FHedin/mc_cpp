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

#include "Angle.hpp"

Angle::Angle() { }

Angle::Angle(int _a1, int _a2, int _a3, int _typ, double _k, double _theta0)
: at1(_a1), at2(_a2), at3(_a3), type(_typ), k(_k), theta0(_theta0) { }

Angle::~Angle() { }

double Angle::getTheta0() const
{
    return theta0;
}

double Angle::getK() const
{
    return k;
}

int Angle::getType() const
{
    return type;
}

int Angle::getAt3() const
{
    return at3;
}

int Angle::getAt2() const
{
    return at2;
}

int Angle::getAt1() const
{
    return at1;
}

std::ostream& operator<<(std::ostream& overloadStream, const Angle& ang)
{
    ang.toString(overloadStream);

    return overloadStream;
}

void Angle::toString(std::ostream& stream) const
{
    stream << std::fixed << std::setprecision(6);
    stream << "Angle" << '\t';
    stream << at1 << '\t' << at2 << '\t' << at3 << type << '\t';
    stream << k << '\t' << theta0;
}

