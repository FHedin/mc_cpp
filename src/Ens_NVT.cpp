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

#include <iostream>

#include "Ens_NVT.h"

Ens_NVT::Ens_NVT(int _N, double _V, double _T)
{
    type = NVT;
    N = _N;
    V = _V;
    T = _T;

    std::cout << "Using ensemble NVT : N = " << N << " V = " << V << " T = " << T << std::endl;
}

Ens_NVT::~Ens_NVT()
{
    //dtor
}

std::string Ens_NVT::whoami()
{
    std::string name("NVT");
    return name;
}

