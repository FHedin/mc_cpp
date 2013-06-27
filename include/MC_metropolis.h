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

#ifndef MC_METROPOLIS_H
#define	MC_METROPOLIS_H

#include "MC.h"

class MC_metropolis : public MC
{
public:
    MC_metropolis(std::vector<Atom>& _at_List, PerConditions& _pbc,
                  Ensemble& _ens, FField& _ff, List_Moves& _mvlist, int _steps, double _dmax,
                  int _update_frequency);
    virtual ~MC_metropolis();

    virtual void run();

private:

    virtual void apply_criterion(int natom, int nmvtyp, int imvtyp, int imvatm);
    
    //    virtual void apply_criterion(Atom const& oldAt, Atom const& newAt, int candidate); //one atom
    //    virtual void apply_criterion(std::vector<Atom>& candidateVector); //all atoms

};

#endif	// MC_METROPOLIS_H

