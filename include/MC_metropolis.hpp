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

//#include "Global_include.hpp"

#include "MC.hpp"

class MC_metropolis : public MC
{
public:
	MC_metropolis(AtomList& _at_List, PerConditions& _pbc, Ensemble& _ens, FField& _ff, List_Moves& _mvlist, int _steps, int _save_freq, uint64_t _seed = 0);
    ~MC_metropolis();

    virtual void run();

private:

    virtual void apply_criterion(double de);

};

#endif	// MC_METROPOLIS_H

