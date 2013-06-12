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

#include "FField.h"

#ifndef FFIELD_MDBAS_H
#define	FFIELD_MDBAS_H

class FField_MDBAS : public FField
{
public:
    FField_MDBAS(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens);
    virtual ~FField_MDBAS();

    // The Lennard-Jones potential
    virtual double getLJ(bool dV); //all atoms
    virtual double getLJ(std::vector<Atom>& candidateVec, bool dV); //all atoms
    virtual double getLJ(Atom const& newAt, int candidate, bool dV); //one atom

private:

};

#endif	/* FFIELD_MDBAS_H */

