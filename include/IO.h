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

#ifndef IO_H
#define	IO_H

#include "Atom.h"
#include "Ensemble.h"
#include "PerConditions.h"

class IO
{
public:
    IO(std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens);
    virtual ~IO();

protected:
    std::vector<Atom>& at_List;
    PerConditions& pbc;
    Ensemble& ens;

    virtual void read_coord() = 0;
    virtual void read_ff() = 0;

};

#endif	/* IO_H */

