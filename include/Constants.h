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

#ifndef CONSTANTS_H
#define CONSTANTS_H

#include "Global_include.h"

namespace CONSTANTS
{
    const double elemchg = 1.602176565e-19;
    const double angstr = 1.e-10;
    const double calory = 4.184;
    const double kcaltoiu = 418.4;
    const double clight = 299792458.;
    const double NA = 6.02214129e+23;
    const double bartoiu = 6.02214129e-3;
    const double kboltz = 1.3806488e-23;
    const double rboltz = 8.3144621;
    const double rboltzui = 0.83144621;
    const double mu0 = 1.e-7;
    const double chgcharmm = 332.0716;
    const double chgnamd = 332.0636;
    const double chgdlpolyiu = 138935.4835;
    const double sq6rt2 = 1.122462048309373;
    const double PI = 3.141592653589793;
    const double TWOPI = 6.283185307179586;
    const double SQRTPI = 1.772453850905516;
    const double watercomp = 0.007372;
}

#endif // CONSTANTS_H