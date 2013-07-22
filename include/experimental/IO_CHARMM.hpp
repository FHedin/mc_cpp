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

#ifndef IO_CHARMM_H
#define	IO_CHARMM_H

#include "Global_include.hpp"

#ifdef CHARMM_EXPERIMENTAL

#include "IO.hpp"

class IO_CHARMM : public IO
{
public:
    IO_CHARMM(std::string configf_name, std::vector<Atom>& _at_List, PerConditions& _pbc, Ensemble& _ens);
    virtual ~IO_CHARMM();

private:
    virtual void read_coord();
    virtual void read_ff();

    FILE* conff;
};

#endif  //CHARMM_EXPERIMENTAL

#endif	/* IO_CHARMM_H */

