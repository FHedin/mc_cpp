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

#ifdef NPT_EXPERIMENTAL

#include <iostream>

#include "Ens_NPT.h"

Ens_NPT::Ens_NPT(int _N, double _P, double _T)
{
    type = NPT;
    N = _N;
    P = _P;
    T = _T;
    
    std::cout << "Using ensemble NPT : N = " << N << " P = " << P << " T = " << T << std::endl;
}

Ens_NPT::~Ens_NPT()
{
}

std::string Ens_NPT::whoami()
{
    std::string name("NPT");
    return name;
}

#endif //NPT_EXPERIMENTAL