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

#include "FField_MDBAS.h"

FField_MDBAS::FField_MDBAS(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens) : FField(_at_List, _pbc, _ens)
{
}

FField_MDBAS::~FField_MDBAS()
{
}

double FField_MDBAS::getLJ(bool dV)
{
    return 0.0;
}

double FField_MDBAS::getLJ(std::vector<Atom>& candidateVec, bool dV)
{
    return 0.0;
}

double FField_MDBAS::getLJ(Atom const& newAt, int candidate, bool dV)
{
    return 0.0;
}
